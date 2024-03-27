from __future__ import print_function
import platform
import sys
import json
import pprint
from collections.abc import Iterable

sys.path.insert(0,'../pycparser')
sys.path.insert(0,'../pycparser/pycparser/')
sys.path.insert(0,'../pycparser/examples/')
sys.path.insert(0,'../pycparser/utils/')
# This is not required if you've installed pycparser into
# your site-packages/ with setup.py
#
sys.path.extend(['.', '..', '../pycparser/pycparser'])

from pycparser import c_ast, parse_file
from IPython import embed

FUNC_CREATE_TASK = "pthread_create"
READ_SHARED_VAR  = "iAutoSyncRead"
WRITE_SHARED_VAR = "iAutoSyncWrite"
READ_TO_UPDATE_SHARED_VAR = "iAutoSyncReadToUpdate"
UPDATE_SHARED_VAR = "iAutoSyncUpdate"
PROCEED_ON_EVENT = "iAutoSyncProceedOnEvent"
                     
PATH_JSON = "../05_Workspace/parser_out.json"

shared_var_usage = {}
auto_sync_calls = {}
intentions = {}           # ToDo: this should be called shared_var_dependencies
general_intentions = {}   # ToDo: this should be called intentions

########################################################################
##                       HELPER FUNCTIONS                             ##
########################################################################
# From StackExchange
def flatten(lis: list) -> list:
     for item in lis:
         if isinstance(item, Iterable) and not isinstance(item, str):
             for x in flatten(item):
                 yield x
         else:
             yield item

def del_duplicates(lis: list) -> list:
    '''
    Delete duplicates of a list
    '''
    return list(dict.fromkeys(lis))


def get_shared_var_from_auto_sync_call(node: c_ast.Node) -> str:
    # arg_pos is the position of the shared-variable in the AutoSync function signature
    if node.name.name == READ_SHARED_VAR or \
       node.name.name == READ_TO_UPDATE_SHARED_VAR:   
       arg_pos = 1 
    elif node.name.name == WRITE_SHARED_VAR or \
        node.name.name == UPDATE_SHARED_VAR:
        arg_pos = 0
    else:
        return ""
        print(f'[PARSER ERROR] AutoSync interface is unknown: {node}')  
        embed()
        exit(1)
        
    try:         
        if isinstance(node.args.exprs[arg_pos], c_ast.ID):
            # Handle cases where the pointer is passed directly (without &)
            shared_var = node.args.exprs[arg_pos].name
        elif isinstance(node.args.exprs[arg_pos], c_ast.StructRef):
            # Handle cases for pointer to structs (StructA->fieldB)
            shared_var = node.args.exprs[arg_pos].name.name + \
                        node.args.exprs[arg_pos].type + \
                        node.args.exprs[arg_pos].field.name
        elif isinstance(node.args.exprs[arg_pos].expr, c_ast.ID):
            # Handle variables with &
            shared_var = node.args.exprs[arg_pos].expr.name
        elif isinstance(node.args.exprs[arg_pos].expr, c_ast.ArrayRef):
                # Handle cases for pointer to structs with arrays (Global->transtimes[0])
            shared_var = node.args.exprs[arg_pos].expr.name.name.name + \
                            node.args.exprs[arg_pos].expr.name.type + \
                            node.args.exprs[arg_pos].expr.name.field.name
        elif isinstance(node.args.exprs[arg_pos].expr, c_ast.StructRef):
            shared_var = node.args.exprs[arg_pos].expr.name.name + \
                        node.args.exprs[arg_pos].expr.type + \
                        node.args.exprs[arg_pos].expr.field.name
        else:
            print(f'[PARSE ERROR] Shared variable of iAutoSyncRead* could not be recognized: {node.args.exprs[arg_pos]}')
            exit(1)

        return(shared_var) 
    except Exception as ex:
        print(f'[PARSER ERROR] Unexpected exception: {ex}')
        print(f'[PARSER ERROR] When parsing this node: {node}')
        exit(1)        


# Get all the existing threads 
class ThreadCreationVisitor(c_ast.NodeVisitor):
    def __init__(self):        
        # main should always exist and is not created with pthread_create
        self.existing_threads = ['main']

    def visit_FuncCall(self, node):
        if node.name.name == FUNC_CREATE_TASK:
            thread = node.args.exprs[2].expr.name
            self.existing_threads.append(thread)        

        # Visit args in case they contain more func calls.
        if node.args:
            self.visit(node.args)

    def get_existing_threads(self):
        return self.existing_threads

    def show(self):
        print(self.existing_threads)


# Get the quantity of each existing thread
class NoOfThreadsVisitor(c_ast.NodeVisitor):
    def __init__(self, existing_threads):
        self.existing_threads = existing_threads
        self.no_of_threads = {}
        self.no_of_threads_outside_loop()


    def no_of_threads_outside_loop(self):
        for thread in self.existing_threads:
            self.no_of_threads[thread] = self.existing_threads.count(thread)


    def visit_For(self, node):
        if FUNC_CREATE_TASK in node.stmt.show.__str__():
            for thread in self.existing_threads:
                if thread in node.stmt.show.__str__():
                    self.no_of_threads[thread] += 1  
    

    def visit_While(self, node):
        if FUNC_CREATE_TASK in node.stmt.show.__str__():
            for thread in self.existing_threads:
                if thread in node.stmt.show.__str__():
                    self.no_of_threads[thread] += 1

    def get_no_of_threads(self):      
        return self.no_of_threads


    def show(self):
        print(self.no_of_threads)


# Get the existing shared-variables 
class SharedVarVisitor(c_ast.NodeVisitor):
    def __init__(self):        
        self.existing_shared_var = []


    def visit_FuncCall(self, node):        
        self.existing_shared_var.append(get_shared_var_from_auto_sync_call(node))

        # Visit args in case they contain more func calls.
        if node.args:
            self.visit(node.args)


    def get_existing_shared_var(self):
        return self.existing_shared_var


    def show(self):
        print(self.existing_shared_var)

# Get coupled variables
class IntentionsVisitor(c_ast.NodeVisitor):
    def __init__(self, intention_var):       
        self.intention_var = intention_var
        self.depends_on = []
        self.constant_init_by_main = []


    def visit_Decl(self, node):
        for var in self.intention_var:  
                
            if var == node.name:
                if "xIntentionRootN" in self.intention_var:
                    pass
                    #embed()
                if node.init is None:
                    print(f"!!! [PARSER INFO] No intention has been specified for {node.name}")
                else:
                    # Iterate over all fields in the xAutoSyncIntentions struct 
                    for intention in node.init.exprs:
                        if intention.name[0].name == "pvDependsOn":                            
                            if isinstance(intention.expr.expr, c_ast.StructRef):
                                dependent_var = intention.expr.expr.name.name + \
                                                intention.expr.expr.type + \
                                                intention.expr.expr.field.name
                            elif isinstance(intention.expr.expr, c_ast.ID):
                                dependent_var = intention.expr.expr.name
                            else:
                                print(f'[PARSE ERROR] Could not read shared-variable dependency: {intention}')
                                exit(1)
                            
                            self.depends_on.append(dependent_var)
                        elif intention.name[0].name == "bConstantInitByMain":
                            if intention.expr.value == '1':
                                # Flag is set to true
                                self.constant_init_by_main.append("bConstantInitByMain") 
                    

    def get_dependecies(self):
        self.depends_on = del_duplicates(self.depends_on)
        return self.depends_on


    def get_constant_init_by_main(self):
        self.constant_init_by_main = del_duplicates(self.constant_init_by_main)
        return self.constant_init_by_main



# Get the existing local and global variables and their types
class VarDeclVisitor(c_ast.NodeVisitor):
    def __init__(self):       
        self.existing_var = {}


    def visit_Decl(self, node):        
        if isinstance(node.type, c_ast.TypeDecl):
            self.existing_var[node.name] = ' '.join(node.type.type.names)
        elif isinstance(node.type, c_ast.ArrayDecl):
            arr_name = node.name
            if isinstance(node.type.type.type, c_ast.IdentifierType):
                arr_type = node.type.type.type.names                
            elif isinstance(node.type.type.type, c_ast.TypeDecl):
                arr_type = node.type.type.type.type.names
            arr_size = node.type.dim.value
            self.existing_var[f'{arr_name}[{arr_size}]'] = ' '.join(arr_type)


    def get_existing_var(self):
        return self.existing_var


    def show(self):
        print(self.existing_var)


# Get the existing shared-variables 
class SharedVarUsageVisitor(c_ast.NodeVisitor):
    def __init__(self, existing_threads):        
        # main should always exist and is not created with pthread_create
        self.shared_var_usage = []
        self.existing_threads = existing_threads


    def visit_FuncDef(self, node):
        # Function Name: node.decl.name
        for thread in existing_threads:
            if node.decl.name == thread:
                v = FuncCallVisitor(thread)
                v.visit(node)  


class FuncCallVisitor(c_ast.NodeVisitor):
    def __init__(self, thread):        
        self.callees = []  
        self.thread = thread
        shared_var_usage[self.thread] = {"Read": list(), 
                                         "Write": list(),
                                         "ReadToUpdate": list(),
                                         "Update": list(),
                                         "Quantity": 0}        


    def visit_FuncCall(self, node):
        self.callees.append(node.name.name)
        line_no = int(node.coord.line) 
        func = node.name.name

        if func == READ_SHARED_VAR:
            #if isinstance(node.args.exprs[1].expr, c_ast.StructRef):
            #    shared_var = node.args.exprs[1].expr.name.name + \
            #                 node.args.exprs[1].expr.type + \
            #                 node.args.exprs[1].expr.field.name
            #elif isinstance(node.args.exprs[1], c_ast.ID):
                # Handle cases where the pointer is passed directly (without &)
            #    shared_var = node.args.exprs[1].name
            #else:
            #    shared_var = node.args.exprs[1].expr.name            
            shared_var = get_shared_var_from_auto_sync_call(node)
            shared_var_usage[self.thread]["Read"].append(shared_var)
            auto_sync_calls[line_no] = (READ_SHARED_VAR, shared_var)
        
            if shared_var not in intentions:
                intentions[shared_var] = []
            intentions[shared_var].append(node.args.exprs[3].name)

        if func == READ_TO_UPDATE_SHARED_VAR:         
            #if isinstance(node.args.exprs[1], c_ast.ID):
                # Handle cases where the pointer is passed directly (without &)
            #    shared_var = node.args.exprs[1].name   
            #elif isinstance(node.args.exprs[1].expr, c_ast.StructRef):
            #    shared_var = node.args.exprs[1].expr.name.name + \
            #                 node.args.exprs[1].expr.type + \
            #                 node.args.exprs[1].expr.field.name            
            #else:
            #    shared_var = node.args.exprs[1].expr.name            
            shared_var = get_shared_var_from_auto_sync_call(node)
            shared_var_usage[self.thread]["ReadToUpdate"].append(shared_var)
            auto_sync_calls[line_no] = (READ_TO_UPDATE_SHARED_VAR, shared_var)
            
            if shared_var not in intentions:
                intentions[shared_var] = []
            intentions[shared_var].append(node.args.exprs[3].name)     
             
            
        if func == WRITE_SHARED_VAR:
            #if isinstance(node.args.exprs[0], c_ast.ID):
                # Handle cases where the pointer is passed directly (without &)
            #    shared_var = node.args.exprs[0].name
            #elif isinstance(node.args.exprs[0].expr, c_ast.StructRef):
            #    shared_var = node.args.exprs[0].expr.name.name + \
            #                 node.args.exprs[0].expr.type + \
            #                 node.args.exprs[0].expr.field.name
            #else:
            #    shared_var = node.args.exprs[0].expr.name            
            shared_var = get_shared_var_from_auto_sync_call(node)
            shared_var_usage[self.thread]["Write"].append(shared_var)
            auto_sync_calls[line_no] = (WRITE_SHARED_VAR, shared_var)

            if shared_var not in intentions:
                intentions[shared_var] = []
            intentions[shared_var].append(node.args.exprs[3].name)  

        if func == UPDATE_SHARED_VAR:
            #if isinstance(node.args.exprs[0], c_ast.ID):
                # Handle cases where the pointer is passed directly (without &)
            #    shared_var = node.args.exprs[0].name
            #elif isinstance(node.args.exprs[0].expr, c_ast.StructRef):
            #    shared_var = node.args.exprs[0].expr.name.name + \
            #                 node.args.exprs[0].expr.type + \
            #                 node.args.exprs[0].expr.field.name
            #else:
            #    shared_var = node.args.exprs[0].expr.name            
            shared_var = get_shared_var_from_auto_sync_call(node)
            shared_var_usage[self.thread]["Update"].append(shared_var)
            auto_sync_calls[line_no] = (UPDATE_SHARED_VAR, shared_var)  

            if shared_var not in intentions:
                intentions[shared_var] = []
            intentions[shared_var].append(node.args.exprs[3].name)        

        if func == PROCEED_ON_EVENT:
           
            event = node.args.exprs[0].name
            no_of_threads = node.args.exprs[1].name
            auto_sync_calls[line_no] = (PROCEED_ON_EVENT, event, no_of_threads)  
        
        # Visit args in case they contain more func calls.
        if node.args:
            self.visit(node.args)

    def get_auto_sync_read_usage(self):
        return self.auto_sync_read 

def show_func_calls(ast, funcname):    
    v = FuncCallVisitor(funcname)
    v.visit(ast)  
    
    
def get_existing_threads(ast) -> list:    
    v = ThreadCreationVisitor()
    v.visit(ast)
    return v.get_existing_threads()

def get_no_of_threads(ast, existing_threads) -> dict:
    v = NoOfThreadsVisitor(existing_threads)
    v.visit(ast)
    return v.no_of_threads


def get_existing_shared_var(ast) -> list:    
    v = SharedVarVisitor()
    v.visit(ast)
    return set(v.get_existing_shared_var())


def get_shared_var_usage(ast, existing_threads) -> list:    
    v = SharedVarUsageVisitor(existing_threads)
    v.visit(ast)
    #return set(v.get_existing_shared_var())

# Why do we need this?
# The FuncCallVisitor is only able to detect the calls "in the first level"
# If we have an AutoSync call inside another function, this is not detected.
# That's why we need this FuncDefVisitor
class FuncDefVisitor(c_ast.NodeVisitor):
    def __init__(self):    
        pass    

    def visit_FuncDef(self, node):
        try:
            line_no = int(node.coord.line) 
            func = node.decl.name

            if func == PROCEED_ON_EVENT:
                event = node.args.exprs[0].name
                no_of_threads = node.args.exprs[1].name
                auto_sync_calls[line_no] = (PROCEED_ON_EVENT, event, no_of_threads)

            if func == UPDATE_SHARED_VAR:  
                embed()       
                shared_var = get_shared_var_from_auto_sync_call(node)
                shared_var_usage[self.thread]["Update"].append(shared_var)
                auto_sync_calls[line_no] = (UPDATE_SHARED_VAR, shared_var)  

                if shared_var not in intentions:
                    intentions[shared_var] = []
                intentions[shared_var].append(node.args.exprs[3].name)         
            
            v = FuncCallVisitor(node.decl.name)
            v.visit(node)
        except Exception as ex:
            print(f'Exception when visiting Function Definitions: {ex}')
            embed()
      



if __name__ == "__main__":
    filename  = sys.argv[1]   
    print(f'MATHEUS: {filename}')
    ast = parse_file(filename, use_cpp=True,
                               cpp_path='gcc',
                               cpp_args=['-E', r'-Iutils/fake_libc_include'])
    
    g = VarDeclVisitor()
    g.visit(ast)
    existing_var = g.get_existing_var()
    #g.show()

    t = FuncDefVisitor()
    t.visit(ast)
            
    # Print information obtained with static analysis
    existing_threads = get_existing_threads(ast)
    no_of_threads = get_no_of_threads(ast, existing_threads)
    existing_shared_var = get_existing_shared_var(ast)
    
    get_shared_var_usage(ast, existing_threads)

    #print(f'Existing threads: {set(existing_threads)}\n')
    #print(f'Quantity of each thread: {no_of_threads}\n')
    #print(f'Existing shared-variables: {set(existing_shared_var)}\n')
    #print(f'Usage of shared-variables: {shared_var_usage}\n')

    for thread in existing_threads:
        shared_var_usage[thread]["Quantity"] = no_of_threads[thread]

    existing_shared_var = []
    for thread, info in shared_var_usage.items():
        for key in info:
            if key == 'Read' or key == 'Write' or \
               key == 'ReadToUpdate' or key == 'Update':
                existing_shared_var.append(info[key])

    if [] in existing_shared_var: existing_shared_var.remove([])
    existing_shared_var = set(list(flatten(existing_shared_var)))

    # Check intentions for plausibility (i.e. check conflicting intentions)
    for key, value in intentions.items():
        intentions_plausible = value.count(value[0]) == len(value)

        if not intentions_plausible:
            print(50*"-")
            print(f"Shared-variable, Intentions: {key, value}")            
            print(50*"-")
            print("        DEBUG INFORMATION     ")
            print(f'[PARSER ERROR] Shared-variable {key} has conflicting intentions!')
            exit(1)        
    
        v = IntentionsVisitor(value)
        v.visit(ast)
        intentions[key] = v.get_dependecies()
        general_intentions[key] = v.get_constant_init_by_main()
        
    print(50*"-")
    parser_output = []
    parser_output.append(shared_var_usage)
    parser_output.append(existing_var)
    parser_output.append(auto_sync_calls)
    parser_output.append(intentions)
    parser_output.append(general_intentions)    

    json_file = json.dumps(parser_output, sort_keys=True, indent=2)
    print(json_file)
    
    with open(PATH_JSON, "w") as output:
        output.write(json_file)

    print(f'-----> Quantity of AutoSync calls: {len(auto_sync_calls.keys())}')

            
    
