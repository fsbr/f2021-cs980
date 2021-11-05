# i need to be able to make a lot of different environments for my thing to open up. 
"""
output: 
int:    10                      # size x
int:    10                      # size y
string: - - - - - - - - - -     # "#" means an obstacle 
string: - # - - - - - - - -     # "-" means a free spacd
string: - - # - - # # - - -     # there are no spaces in the real input map
string: - - - - - - # - - -
string: - - - - - - - - - -
string: - - - - - - - - - -
string: - - - - - - - - - -
string: - - - - - - - - - -
string: - # - - - # # - - -
string: - - - - - - - - - -
float:  start x
float:  start y
float:  goal x
float:  goal y
"""
import math
import random
def make_environment(sizeX=10, sizeY=10, startX=0.25, startY=0.25, goalX=9.5, goalY=9.5, environmentIdx = 0):
    # indexing thing
    sizeXidx = sizeX - 1
    sizeYidx = sizeY - 1
    f = open("grid_envs/environment%s.txt"%environmentIdx, "w")
    # 
    f. write(str(sizeX) + "\n")
    f. write(str(sizeY) + "\n")
    # first step, assemble a string of ten characters 
    rowString = []
    for j in range(sizeY):
        for i in range(sizeX):
            #print(" i is ", i)
            randomNumber = random.uniform(0,1)
            #print(randomNumber) 
            if (j == math.floor(sizeY-startY)) and (i == math.floor(startX)):  
                rowString.append("-")
            elif (j== math.floor(sizeY - goalY)) and (i== math.floor(goalX)): #start y
                rowString.append("-")
            elif (randomNumber < 0.20):
                rowString.append("#")
            else:
                rowString.append("-")
        #print(rowString)
        # instead of printing the rowstring, we want to write it
        for character in rowString:
            f.write(character)
            f.write(" ")
        f.write("\n")
        rowString = []

    f.write(str(startX) + "\n")
    f.write(str(startY) + "\n")
    f.write(str(goalX) + "\n")
    f.write(str(goalY) + "\n")
    f.write("\n")
    f.close()
if __name__ == "__main__":
    random.seed(31415)
    sizeX = 10
    sizeY = 10
    for i in range(10000):
        startX = random.uniform(0,sizeX)
        startY = random.uniform(0,sizeY)
        goalX = random.uniform(0,sizeX)
        goalY = random.uniform(0,sizeY)
        make_environment(sizeX, sizeY, startX, startY, goalX, goalY, environmentIdx=i)
    #make_environment()
    
