import os
newfile = ''
for object in os.listdir():
    if object.find("GEMCSCLUT")==-1:
        continue
    with open(object) as file:
        lines = file.readlines()
        for line in lines:
            try:
                superchambereta = int(line.split(" ")[0])
                chamber = superchambereta//10
                eta = superchambereta%10
                newkey = 8*(chamber-1)+(eta-1)
                line = line.replace(str(superchambereta),str(newkey))
            except:
                pass
            newfile += line
    with open(object, "w") as outfile:
        outfile.write(newfile)