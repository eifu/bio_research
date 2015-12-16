## bio_proj3_1.py
##
## by eifu -- Oct.29
## mentored by Dr. Regan
##
##
##
"""
this file is to create pdf files of T status, the average value of when the status get turned from OFF to ON.
"""
import matplotlib.pyplot as pyplot
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import histone

##############################

NUM_OF_HISTONE = 81
BEFORE_PROMOTER = 40
WIDTH = 10
TIME = 500

##############################
def createRandomHistoneList(percentage):
    dstList = []
    for i in range(NUM_OF_HISTONE):
        if(i==0):dstList.append(histone.Histone(i-BEFORE_PROMOTER,percentage=percentage))
        else:
            dstList.append(histone.Histone(i-BEFORE_PROMOTER,preNode=dstList[i-1],percentage=percentage))
            dstList[i-1].set_nextNode(dstList[i])

    return dstList

def trackingHistone(histoneList,R,A,T,Eext):
    trackerList = [[] for i in range(NUM_OF_HISTONE)]
    """
    tracker List is a two dimentional list,
    intuitively, the first element in the trackerList is a list that only tracks the 1st histone chronically, and
    the second element in the trackerList is a list that only tracks the 2nd histone.
    """
    TEextTrackerList=[]
    TEextTrackerList.append((T,Eext))
    start = BEFORE_PROMOTER - WIDTH/2
    end = BEFORE_PROMOTER + WIDTH/2
    check = 0
    for t in range(TIME):
        for i in range(NUM_OF_HISTONE):
            trackerList[i].append(histoneList[i].condition)
            if(start <= i and i<= end):
                check += histoneList[i].condition
            histoneList[i].k_plus()
            histoneList[i].k_minus()
        T = A and (check == 0)
        Eext = ((not T) and (not A)) or R
        TEextTrackerList.append((T,Eext))

        check = 0
        if(Eext == True):
            histoneList[BEFORE_PROMOTER].condition = 1

    return trackerList,TEextTrackerList


def main():
    progress = 0
    R = 0
    A = 1
    for T in range(2):
        for Eext in range(2):
            for percentage in range(0,101,50):
                submain(R,A,T,Eext,percentage)

                progress += 25/3
                print ("{} % done".format(progress))

    # submain(0,0,1,1,50)

def takeListAndConvert(ListForResult,ListtoBeAdded):
    result = [ListForResult[i] + ListtoBeAdded[i] for i in range(TIME)]
    return result

def submain(R,A,T,Eext,percentage):


    time = np.linspace(0,TIME-1,TIME)
    fig = pyplot.figure()
    title = "exp4/exp4_R{}_A{}_T{}_Eext{}_percentage{}_with200histones.pdf".format(R,A,T,Eext,percentage)
    pp = PdfPages(title)

    histoneList = createRandomHistoneList(percentage)
    trackerList,TEextTrackerList = trackingHistone(histoneList,R,A,T,Eext)
    T_List = [TEextTrackerList[i][0] for i in range(TIME)]

    for i in range(199):
        histoneList = createRandomHistoneList(percentage)
        trackerList,TEextTrackerList = trackingHistone(histoneList,R,A,T,Eext)
        T_toBeAdded_List = [TEextTrackerList[i][0] for i in range(TIME)]
        T_List = takeListAndConvert(T_List,T_toBeAdded_List)

    submain2(fig,time,T_List,R,A,T,Eext,percentage)
    pp.savefig(fig)
    pp.close()
    #pyplot.show()


## this function is to create a graph of T status.
## T is one of TFs, that enable
def submain2(fig,time,T_List,R,A,T,Eext,percentage):
    cx = fig.add_subplot(111)
    # T_List = T_list_maker(trackerList)
    #T_List = [TEextTrackerList[i][0] for i in range(TIME)]
    cx.plot(time,T_List,"-",color="red")
    #cx.scatter(time,T_List)
    cx.set_xlabel("time")
    cx.set_xlim(-0.5,TIME)
    cx.set_ylabel("T")
    cx.set_yticks((0,200))
    cx.set_yticklabels(("Off","On"))
    cx.set_ylim(-5,205)

    cx.set_title("R={}, A={}, T={}, Eext={}, percentage={}".format(R,A,T,Eext,percentage))



if __name__ == "__main__":
    main()
