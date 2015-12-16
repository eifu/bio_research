## bio_proj3_3.py
##
## eifu by Nov. 19 -- mentored by Dr. Regan
##
"""
this file is to create a pdf files with methylated, unmethylated, and acetilated histones.
"""


import histone
from numpy.random import sample
import matplotlib.pyplot as pyplot
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
##
NUM_OF_HISTONE = 81
BEFORE_PROMOTER = 40
TIME = 500
WIDTH = 10
##


def main():

    progress = 0
    for R in range(2):
        for A in range(2):
            for T in range(2):
                for Eext in range(2):

                    submain(R,A,T,Eext,0)
                    submain(R,A,T,Eext,50)
                    submain(R,A,T,Eext,100)
                    progress += 6.25
                    print("{} %% done".format(progress))

    # submain(0,0,1,1,50)

def submain(R,A,T,Eext,percentage):

    histoneList = histone.createRandomHistoneList(percentage,A)
    trackerList, TEextTrackerList = histone.trackingHistone(histoneList,R,A,T,Eext,TIME=TIME)


    time = np.linspace(0,TIME-1,TIME)
    fig = pyplot.figure()
    title = "exp7.2/exp7_R{}_A{}_T{}_Eext{}_percentage{}.pdf".format(R,A,T,Eext,percentage)
    pp = PdfPages(title)



    submain1(fig,time,TEextTrackerList)
    submain2(fig,time,TEextTrackerList)
    submain3(fig,trackerList,R,A,T,Eext,percentage)

    pp.savefig(fig)
    pp.close()
    # pyplot.show()



def submain1(fig,time,TEextTrackerList):
    """
    this function is to create a graph of Eext status
    """

    ax = fig.add_subplot(514)
    # print(TIME)
    # print(len(TEextTrackerList))
    EextList = [TEextTrackerList[i][1] for i in range(TIME)]
    ax.plot(time,EextList,"-",color="blue",drawstyle="steps")
    ax.tick_params(bottom = "off",labelbottom='off')
    #ax.scatter(time,trackerList[BEFORE_PROMOTER])
    ax.set_xlim(-0.5,TIME)
    #ax.set_xticks([i for i in range(0, TIME+1) if (i%5==0) ])
    ax.set_ylabel("Eext")
    ax.set_yticks((0,1))
    ax.set_yticklabels(("Off","On"))
    ax.set_ylim(-0.5,1.5)


def submain2(fig,time,TEextTrackerList):
    """
    this function is to create a graph of T status
    """
    cx = fig.add_subplot(515)
    # T_List = T_list_maker(trackerList)
    T_List = [TEextTrackerList[i][0] for i in range(TIME)]
    cx.plot(time,T_List,"-",color="red",drawstyle="steps")
    #cx.scatter(time,T_List)
    cx.set_xlabel("time")
    cx.set_xlim(-0.5,TIME)
    cx.set_ylabel("T")
    cx.set_yticks((0,1))
    cx.set_yticklabels(("Off","On"))
    cx.set_ylim(-0.5,1.5)


def submain3(fig,trackerList,R,A,T,Eext,percentage):
    """
    this function is to create a graph of all histones status
    """
    bx = fig.add_subplot(211)
    bx.tick_params(left ="off",labelleft="off")

    for h in range(NUM_OF_HISTONE):
        trackerM = [i for i in range(TIME) if trackerList[h][i] == "m"]
        y = [NUM_OF_HISTONE-1-h for i in range(len(trackerM))]
        bx.plot(trackerM,y,",",color = "blue")
        trackerA = [i for i in range(TIME) if trackerList[h][i] == "a"]
        y = [NUM_OF_HISTONE-1-h for i in range(len(trackerA))]
        bx.plot(trackerA,y,",",color = "red")



    bx.set_xlim(-0.5,TIME)
    bx.set_ylabel("histones' status")

    bx.set_ylim(-1,NUM_OF_HISTONE+0.5)
    bx.set_title("R={}, A={}, T={}, Eext={}, percentage={}".format(R,A,T,Eext,percentage))



if __name__ == "__main__":
    main()
