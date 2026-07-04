import os
newfile = ''
for object in os.listdir("../AlignmentCorrection_0/"):
    if object.find("GEMCSCLUT")==-1:
        continue
    if object.find("4d")==-1:
        continue
    twod = {}
    with open(object.replace("_4d",""), "r") as file:
        lines = file.readlines()
        for line in lines:
            twod[line.split(" ")[0]] = line.split(" ")[1].replace("\n","")
    if object.find("negative")!=-1:
        endcap = "negative"
    else:
        endcap = "positive"
    L1 = {}
    with open("additional_corr_"+endcap+"_endcap_L1.txt", "r") as file:
        lines = file.readlines()
        for line in lines:
            try:
                superchambereta = int(line.split(" ")[0])
                chamber = superchambereta//10
                eta = superchambereta%10
                newkey = 8*(chamber-1)+(eta-1)
                line = line.replace(str(superchambereta),str(newkey))
                L1[line.split(" ")[0]] = line.split(" ")[1].replace("\n","")
            except:
                pass
    L2 = {}
    with open("additional_corr_"+endcap+"_endcap_L2.txt", "r") as file:
        lines = file.readlines()
        for line in lines:
            try:
                superchambereta = int(line.split(" ")[0])
                chamber = superchambereta//10
                eta = superchambereta%10
                newkey = 8*(chamber-1)+(eta-1)
                line = line.replace(str(superchambereta),str(newkey))
                L2[line.split(" ")[0]] = line.split(" ")[1].replace("\n","")
            except:
                pass
    with open("../AlignmentCorrection_0/"+object,"r") as file:
        if object.find("ME11")!=-1:
            chamberm = 3072
            rollm=384
            layerm=192
        else:
            chamberm = 12288
            rollm=768
            layerm=384
        lines = file.readlines()
        for line in lines:
            try:
                key = int(line.split(" ")[0])
                chamber = key // chamberm + 1
                r = key % chamberm

                roll = r // rollm + 1
                r = r % rollm
                
                layer = r // layerm + 1
                pad = r % layerm
                twodkey = str(8*(chamber-1)+(roll-1))
                newcorrection = int(twod[twodkey])
                if layer == 1:
                    newcorrection += int(L1[twodkey])
                else:
                    newcorrection += int(L2[twodkey])
                newfile += line.split(" ")[0]+" "+str(newcorrection)+'\n'
            except:
                newfile += line
    with open(object, "w") as outfile:
        outfile.write(newfile)
