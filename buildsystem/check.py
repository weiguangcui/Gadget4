import re 
import sys
import os

def parseIf(string, defines, fin):
    s = string
    
    while len(s) > 0:
        if s.startswith("defined"):
            s = s[7:]
            continue
            
        if s.startswith("//"):
            return
            
        m = re.search("[a-zA-Z_][a-zA-Z_0-9]*",s)
        if m is not None and m.start() == 0:
            defines.update([m.group()])
            #print "%s : %s"%(string,m.group())
            s = s[m.span()[1]:]
            continue
        
        if s.startswith("/*"):
            m = re.search("\*/",s)
            if m is not None:
                s = s[m.span()[1]:]
                continue
            else:
                return
            
        if s[0] == "\n":
            return
            
        if s == "\\\n":
            string = fin.readline()
            s = string
            continue
            
        if s[0] in '<>+-*/=!&|() 0123456789\t':
            s = s[1:]
            continue
            
        print("Strange character in '%s' detected: '%s', skipping it."%(string, s[0]))
        s = s[1:]
            
            

#Find all macros used  in a .c or .h files
def filter_code(fin):
    defines = set()
    
    first_encountered = False

    line = fin.readline()
    while line != "":
        s = line.lstrip()

        if s.startswith("#include"):
            if first_encountered == False:
                m = re.search("gadgetconfig",s)
                if m is not None:
                    first_encountered = True
                else:
                    print("First header file included ('%s') is not gadgetconfig.h -- please change this\n"%s.rstrip())
                    exit(1)

        if s.startswith("#if "):
            parseIf(s[4:],defines,fin)
        elif s.startswith("#elseif "):
            parseIf(s[8:],defines,fin)
        elif s.startswith("#elif "):
            parseIf(s[6:],defines,fin)    
            
        elif s.startswith("#ifdef ") or s.startswith("#ifndef "):
            if s.startswith("#ifdef "):
                s = s[7:]
            else:
                s = s[8:]
                
            s = s.lstrip()
            m = re.search("[a-zA-Z_][a-zA-Z_0-9]*",s)
            if m is not None and m.start() == 0:
                defines.update([m.group()])
            else:
                print("Strange #ifdef/#ifndef: '%s'. Skipping it.",s)
            
        line = fin.readline()
    
    return defines
    
#Find all items of Template-Config.sh
def filter_template_config(fin):
    defines = set()
    for line in fin:
        s = line.split()
        if(len(s)>0):
            d = re.findall("^#*([a-zA-Z_][a-zA-Z_0-9]*)",s[0])
            for dd in d:
                defines.update([dd])
                #print dd
    
    return defines

#Find all items of src/data/allvars.cc
def filter_template_ioparam(fin):
    defines = set()
    s = fin.read()
    d = re.findall("add_param\(\"([a-zA-Z_][a-zA-Z_0-9]*)\"",s)
    for dd in d:
        defines.update([dd])

    defines.update(["SofteningComovingClass0"])
    defines.update(["SofteningMaxPhysClass0"])
    defines.update(["SofteningClassOfPartType0"])
        
    return defines


#Find all active items on Config.sh    
def filter_config(fin):
    defines = set()
    for line in fin:
        s = line.split()
        if(len(s)>0):
            d = re.findall("^([a-zA-Z_][a-zA-Z_0-9]*)",s[0])
            for dd in d:
                defines.update([dd])
    
    return defines
    
#Find all macros used in Makefile
def filter_makefile(fin):
    defines = set()
    for line in fin:
        s = line.strip()
        if s.startswith("ifeq"):
            d = re.findall("ifeq\s*\(([a-zA-Z_][a-zA-Z_0-9]*)\s*,\s*\$\(findstring",s)
            for dd in d:
                defines.update([dd])
        if s.startswith("ifneq"):
            d = re.findall("ifneq\s*\(([a-zA-Z_][a-zA-Z_0-9]*)\s*,\s*\$\(findstring",s)
            for dd in d:
                defines.update([dd])
    
    return defines
    
#Parse all documented options from README-Config.md
def filter_readme_config(fin):
    defines = set()
    s = fin.read()
#    d = re.findall(r"^\*\*([a-zA-Z_][a-zA-Z_0-9]*)\*\*\s*\n\s*\n[\W\w]{40}",s)
    d = re.findall(r"\*\*([a-zA-Z_][a-zA-Z_0-9]*)\*\*",s)
    for dd in d:
        defines.update([dd])
        
    return defines
    
    
def load(fin):
    s = set()
    
    for i in fin:
        if not i.startswith("#"):
            s.update([i.strip()])
        
    return s
    
def write(s, fout):
    fout = open(fout,'w')
    for i in sorted(s):
        fout.write(i + os.linesep)
        
    fout.close()
    
#Check source files for illegal macros
def check_code(fin, fout, template, extra):
    allowed = filter_template_config(template)
    #print allowed
    allowed.update(filter_config(extra))
    
    used = filter_code(fin)
    #print used
    
    diff =  sorted(used.difference(allowed))
    
    if len(diff) > 0:
        print("\nIllegal macros/options detected in file %s.\nCheck for potential typos and add them either to the file 'Template-Config.sh' or 'defines_extra'\n('defines_extra' is for macros which are either internally defined or should not appear in Template-Config.sh).\nIn case you want to suppress this check, build with 'make build' instead of 'make'.\n"%fin.name)
        for i in diff:
            print(i)
        print("")
        exit(1)
        
    write(used,fout)
    exit(0)
    
#Check Makefile for illegal options
def check_makefile(fin, fout, template, extra):
    allowed = filter_template_config(template)
    #print allowed
    allowed.update(filter_config(extra))
    
    used = filter_makefile(fin)
    #print used
    
    diff =  sorted(used.difference(allowed))
    
    if len(diff) > 0:
        print("\nIllegal macros/options detected in file %s.\nCheck for potential typos and add them either to the file 'Template-Config.sh' or 'defines_extra'\n('defines_extra' is for macros which are either internally defined or should not appear in Template-Config.sh).\nIn case you want to suppress this check, build with 'make build' instead of 'make'.\n"%fin.name)
        for i in diff:
            print(i)
        print("")
        exit(1)
        
    write(used,fout)
    exit(0)
    
#Check Config.sh for illegal options
def check_config(fin, fout, args, extra):
    allowed = filter_config(extra)
    
    for file in args:
        allowed.update(load(open(file,'r')))

    used = filter_config(fin)
    
    diff =  sorted(used.difference(allowed))
    
    if len(diff) > 0:
        print("\nThe following options are active in %s, but are not used in any of the\nsource code files being compiled into the final executable.\nPlease check for typos and deactivate the options.\nIn case you want to suppress this check, build with 'make build' instead of 'make'.\n"%fin.name)
        for i in diff:
            print(i)
        print("")
        exit(1)
        
    write(used,fout)
    exit(0)

#Check whether all Template-=Config options are documented
def check_documentation(fin, fout, fdoc):
    documented = filter_readme_config(fdoc)
    
    used = filter_template_config(fin)
    
    diff =  sorted(used.difference(documented))
    
    ex = False
    
    if len(diff) > 0:
        print("\nERROR: The following options are undocumented in %s, but appear in %s.\n       Please add a proper documentation!\n"%(fdoc.name,fin.name))
        for i in diff:
            print(i)
        print("")
        ex = True
        
    diff =  sorted(documented.difference(used))
    
    if len(diff) > 0:
        print("\nERROR: The following options are documented in %s, but are not used in %s anymore.\n       Please remove redundant documentation!\n"%(fdoc.name,fin.name))
        for i in diff:
            print(i)
        print("")
        ex = True
        
    if ex:
        exit(1)
        
    write(used,fout)
    exit(0)    

#Check whether all src/io/paramaters.c options are documented
def check_parameters(fin, fout, fdoc):
    documented = filter_readme_config(fdoc)
    
    used = filter_template_ioparam(fin)
    
    diff =  sorted(used.difference(documented))
    
    ex = False
    
    if len(diff) > 0:
        print("\nERROR: The following parameters are undocumented in %s, but appear in %s.\n       Please add a proper documentation!\n"%(fdoc.name,fin.name))
        for i in diff:
            print(i)
        print("")
        ex = True
        
    diff =  sorted(documented.difference(used))
    
    if len(diff) > 0:
        print("\nERROR: The following parameters are documented in %s, but are not used in %s.\n       Please remove redundant documentation!\n"%(fdoc.name,fin.name))
        for i in diff:
            print(i)
        print("")
        ex = True
        
    if ex:
        exit(1)
        
    write(used,fout)
    exit(0)    


if __name__ == "__main__":
    if len(sys.argv) < 3:
        exit(1)
        
    mode = int(sys.argv[1])
    if mode < 1 or mode > 5:
        print("Unknown mode")
        exit(1)
        
    fin = open(sys.argv[2],'r')
    fout = sys.argv[3]
    
    if mode == 1:
        print("Checking %s for illegal define macros"%sys.argv[2])
        template = open(sys.argv[4],'r')
        extra = open(sys.argv[5],'r')
        
        check_code(fin, fout, template, extra)
        
    if mode == 2:
        print("Checking active options of %s"%sys.argv[2])
        extra = open(sys.argv[4],'r')
        check_config(fin, fout, sys.argv[5:], extra)
        
    if mode == 3:
        print("Checking %s for illegal define macros"%sys.argv[2])
        template = open(sys.argv[4],'r')
        extra = open(sys.argv[5],'r')
        
        check_makefile(fin, fout, template, extra)
        
    if mode == 4:
        print("Checking %s for documented options"%sys.argv[2])
        template = open(sys.argv[2],'r')
        doc = open(sys.argv[4],'r')
        
        check_documentation(fin, fout, doc)
                
    if mode == 5:
        print("Checking %s for documented parameters"%sys.argv[2])
        template = open(sys.argv[2],'r')
        doc = open(sys.argv[4],'r')
        
        check_parameters(fin, fout, doc)
          
