## histone.py
##
## by eifu -- Nov. 29
## mentored by Dr. Regan
##
##
##########################

from numpy.random import sample


class Histone(object):
    """
    Histone is a python mudule designed for researching
    methylation modeling. This class allows you to create
    a hitone object.
    """
    K_PLUS_DEFAULT= 0.12
    K_PLUS2_DEFAULT= 0.12
    K_MINUS_DEFAULT= 0.117
    K_ACE_DEFAULT = 0.12

    def __init__(self,position=0,
                 K_PLUS=K_PLUS_DEFAULT,
                 K_PLUS2=K_PLUS2_DEFAULT,
                 K_MINUS=K_MINUS_DEFAULT,
                 A_bool = False,
                 K_ACE=K_ACE_DEFAULT,
                 nextNode=None,preNode=None,
                 percentage=50,copy=False,copy_histone=None):
        """
        Initialize the histone. The position must be
        specified, but K_PLUS, K_MINUS, nextNode, preNode,
        and percentage can be specified. If they are not,
        default values variables will be assigned to those
        instance variables.

        position  ... a position of histone, to differentiate
                  from other histone instances.

        A_bool ... a boolean value, activator of transcription factors
                   if A is true, then the K_ACE is the given paramater.
                   else, K_ACE is 0, which means no histones can be turned to be acetilated.
        preNode ... a neighbor of the instances, judged by position - 1
        nextNode ... a neighbor of the instances, judged by position + 1
        """

        if(copy):
            self.position = copy_histone.position

            self.preNode = copy_histone.preNode  # doubly-linked
            copy_histone.preNode.nextNode = self

            self.nextNode = copy_histone.nextNode # doubly-linked
            copy_histone.nextNode.preNode = self



        else:
            self.position = position
            self.preNode = preNode
            self.nextNode = nextNode
            if(A_bool):Histone.K_ACE = K_ACE
            else:Histone.K_ACE = 0
            Histone.K_PLUS = K_PLUS
            Histone.K_PLUS2 = K_PLUS2
            Histone.K_MINUS = K_MINUS



    def set_adjHistone(self,nextNode):
        """
        Set the input to self.nextNode, one of its instance variables.
        """
        self.nextNode = nextNode
        nextNode.preNode = self

    def k_plus(self):
        """
        the unmethylated histone will get methylated by K_PLUS probability
        return methylated object if the histone will get methylated
        """
        return self

    def k_minus(self):
        """
        the methylated histone will be get unmethylated by K_MINUS probability
        or
        the acetilated histone will be get unacetylated by K_MINUS probability
        return unmethylated histone object if the histone will get unmethylated
        """
        if(sample()<Histone.K_MINUS):return UHistone(copy=True,copy_histone=self)
        return self


    def k_ace(self):
        """
        the unmethylated hitone will be get acetilated by K_ACE probability
        return acetilated histone if the histone will get acetilated
        """
        return self


    def display(self):
        sentence = "address: {}\nposition: {}\nstatus: {}\n".format(self,self.position,self.status)

        try:
            sentence += "Prev Histone:{}\n".format(self.preNode)
        except(AttributeError):
            sentence += "Prev Histone: None\n"
        try:
            sentence += "Next Histone: {}".format(self.nextNode)
        except(AttributeError):
            sentence += "Next Histone: None"
        finally:
            print(sentence)
            print("###\n")


class MHistone(Histone):
    def __init__(self,**kwarg):
        super().__init__(**kwarg)
        self.status="m"


class UHistone(Histone):
    def __init__(self,**kwarg):
        super().__init__(**kwarg)
        self.status="u"

    def k_plus(self):
        if(self.preNode.status == "m" and sample() < Histone.K_PLUS):return MHistone(copy=True,copy_histone=self)
        if(self.nextNode.status == "m" and sample() < Histone.K_PLUS):return MHistone(copy=True,copy_histone=self)
        return self

    def k_minus(self): return self

    def k_ace(self):
        if(sample()<Histone.K_ACE):return AHistone(copy=True,copy_histone=self)
        return self

class AHistone(Histone):
    def __init__(self,**kwarg):
        super().__init__(**kwarg)
        self.status="a"

    def k_plus(self):
        if(self.preNode.status == "m" and sample() < Histone.K_PLUS2):return MHistone(copy=True,copy_histone=self)
        if(self.nextNode.status == "m" and sample() < Histone.K_PLUS2):return MHistone(copy=True,copy_histone=self)
        return self


def createRandomHistoneList(percentage=50,A=1,
                            NUM_OF_HISTONE=81,BEFORE_PROMOTER=40,
                            K_PLUS=0.12,
                            K_MINUS=0.117,
                            K_PLUS2=0.12,
                            K_ACE=0.12):
    """
    percentage ... the probability of having methylated hitone.
    """
    dstList = []
    ratio = percentage/100  ## ratio should be float number between 0 and 1

    for i in range(NUM_OF_HISTONE):

        if(sample() < ratio):dstList.append(MHistone(position=i-BEFORE_PROMOTER,
                                                     K_PLUS=K_PLUS,
                                                     K_PLUS2=K_PLUS2,
                                                     K_MINUS=K_MINUS,
                                                     K_ACE=K_ACE,
                                                     A_bool=A,percentage=percentage))
        else:dstList.append(AHistone(position=i-BEFORE_PROMOTER,
                                     K_PLUS=K_PLUS,
                                     K_PLUS2=K_PLUS2,
                                     K_MINUS=K_MINUS,
                                     K_ACE=K_ACE,
                                     A_bool=A,percentage=percentage))
        # dstList[i-1].display()


        dstList[i-1].set_adjHistone(dstList[i])
    dstList[NUM_OF_HISTONE-1].set_adjHistone(dstList[0]) ## Connect the head to tail
    return dstList


def trackingHistone(histoneList,R,A,T,Eext,NUM_OF_HISTONE=81,BEFORE_PROMOTER=40,WIDTH=10,TIME=100):
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
            trackerList[i].append(histoneList[i].status)
            if(start <= i and i<= end):
                if(histoneList[i].status != "a"):
                    check += 1
            temp = histoneList[i]
            histoneList[i] = histoneList[i].k_minus()
            histoneList[i] = histoneList[i].k_ace()
            histoneList[i] = histoneList[i].k_plus()
            # print((temp.status,histoneList[i].status))
        # print (check)
        T = A and (check == 0)
        Eext = ((not T) and (not A)) or R
        TEextTrackerList.append((T,Eext))

        check = 0
        if(Eext == True):
            histoneList[BEFORE_PROMOTER] = MHistone(copy=True,copy_histone=histoneList[BEFORE_PROMOTER])
    # for i in range(NUM_OF_HISTONE):
    #     print(histoneList[i].display())
    return trackerList,TEextTrackerList
