import numpy as np
import matplotlib.pylab as plt
from math import exp


class IsingModel:
    def __init__(self,L:int,step:int,temp:float) -> None:
        #L是长度，size是L*L,步数是step
        self.size=L*L
        self.length=L
        self.nsteps=step
        self.temp=temp
        #温度范围是0.015, 0.015, 4.5
        #self.temp=np.linspace(0.015,4.5,300)
        #随机初始化晶格[-1,1] 
        self.lattice=np.random.choice([-1,1],size=(L,L))
        self.energy=0.0
        self.energy_square=[]
        self.energy_average=0.0
        self.C=0.0
        self.M_average=0.0
        self.M_per_step=[]
        #self.M=0.0
        self.X=0.0
    def Calc_C(self):
        #计算C
        self.C=(sum(self.energy_square)/len(self.energy_square)-self.energy_average**2)/(self.size*1*self.temp)
        return self.C
    def Calc_M(self):
        self.M=self.M_average
        return self.M
    def Calc_X(self):
        self.X=(sum([M**2 for M in self.M_per_step])/(self.nsteps)-self.M_average**2)/self.size
        return self.X
    #def Init(self):
        #self.lattice=np.random.choice([-1,1],size=(self.length,self.length))
    def CalcEnergy(self):
        #E=-sigma_ij 1*i*j 
        self.energy=0.0
        #由于没有外场，J=1，能量简化很多
        #能量只计算最近邻相互作用,采用周期性边界条件
        for i in range(self.length):
            #0-9 if L=10
            for j in range(self.length):
                if i!=0 and i!=self.length-1 and j !=0 and j!=self.length-1:
                    #除边界外
                    self.energy -=  self.lattice[i,j]*self.lattice[i,j+1]+ \
                                    self.lattice[i,j]+self.lattice[i,j-1]+ \
                                    self.lattice[i,j]*self.lattice[i-1,j]+ \
                                    self.lattice[i,j]*self.lattice[i+1,j]
                #处理边界
                elif i==0 and j != 0 and j != self.length-1:
                    #第一行除左右两点
                    self.energy -= self.lattice[i,j]*self.lattice[i,j+1]+\
                                   self.lattice[i,j]*self.lattice[i,j-1]+\
                                   self.lattice[i,j]*self.lattice[i+1,j]+\
                                   self.lattice[i,j]*self.lattice[self.length-1,j]
                elif i==0 and j == (self.length-1):
                    #右上角顶点
                    self.energy -= self.lattice[i,j]*self.lattice[i-1,j]+\
                                   self.lattice[i,j]*self.lattice[i,0]+\
                                   self.lattice[i,j]*self.lattice[self.length-1,j]+\
                                   self.lattice[i,j]*self.lattice[i+1,j]
                elif i==0 and j==0:
                    #左上角顶点
                    self.energy -= self.lattice[0,0]*self.lattice[0,1]+\
                                    self.lattice[0,0]*self.lattice[0,self.length-1]+\
                                    self.lattice[0,0]*self.lattice[1,0]+\
                                    self.lattice[0,0]*self.lattice[self.length-1,0]
                elif i==self.length-1 and j!=0 and j!=self.length-1:
                    #底边除左右两端
                    self.energy -= self.lattice[i,j]*self.lattice[i,j+1]+\
                                   self.lattice[i,j]*self.lattice[i,j-1]+\
                                   self.lattice[i,j]*self.lattice[i-1,j]+\
                                   self.lattice[i,j]*self.lattice[0,j]
                elif i==self.length-1 and j==0:
                    #左下角
                    self.energy -= self.lattice[i,0]*self.lattice[i,1]+\
                                    self.lattice[i,0]*self.lattice[i,self.length-1]+\
                                    self.lattice[i,0]*self.lattice[i-1,0]+\
                                    self.lattice[i,0]*self.lattice[0,0]
                elif i==self.length-1 and j==self.length-1:
                    self.energy -= self.lattice[i,j]*self.lattice[i,0]+\
                                    self.lattice[i,j]*self.lattice[i,j-1]+\
                                    self.lattice[i,j]*self.lattice[i-1,j]+\
                                    self.lattice[i,j]*self.lattice[0,j]
                elif i != 0 and i!=self.length-1 and j==0:
                    #左边除上下两点
                    self.energy -= self.lattice[i,j]*self.lattice[i,j+1]+\
                                    self.lattice[i,j]*self.lattice[i,self.length-1]+\
                                    self.lattice[i,j]*self.lattice[i-1,j]+\
                                    self.lattice[i,j]*self.lattice[i+1,j]
                elif i != 0 and i!=self.length-1 and j==self.length-1:
                    #右边除上下两点
                    self.energy -= self.lattice[i,j]*self.lattice[i,0]+\
                                    self.lattice[i,j]*self.lattice[i,j-1]+\
                                    self.lattice[i,j]*self.lattice[i-1,j]+\
                                    self.lattice[i,j]*self.lattice[i+1,j]    
        return self.energy/4
        #重复计数，是4?  
    #下面是chatGPT给出的能量计算代码
    # def CalcEnergy(self):
    #     self.energy = 0.0
    #     for i in range(self.length):
    #         for j in range(self.length):
    #             nn_sum = self.lattice[(i + 1) % self.length, j] + self.lattice[(i - 1) % self.length, j] \
    #                  + self.lattice[i, (j + 1) % self.length] + self.lattice[i, (j - 1) % self.length]
    #             self.energy -= self.lattice[i, j] * nn_sum
    #     return self.energy / 2
         
    def DrawModel(self):
        fig, ax = plt.subplots()
        cmap = plt.cm.get_cmap('gray')# 'viridis', 'magma', 'jet', 'cool', 'hot', 'spring', 'summer', 'autumn', 'winter'
        im = ax.imshow(self.lattice, cmap=cmap)
        #cbar = ax.figure.colorbar(im, ax=ax)
        ax.set_title('{}*{} Ising model'.format(self.length,self.length))
        ax.set_xticks(np.arange(-0.5, self.length, 1), minor=True)
        ax.set_yticks(np.arange(-0.5, self.length, 1), minor=True)
        ax.grid(which='minor', color='black', linestyle='-', linewidth=1)
        plt.show()
    def MCmove(self):
        k=1
        # k_B， Boltzmann const
        acc_probability=0.0
        step=1
        energy_sum=0.0
        #energy_square=[]
        while step<= self.nsteps:
            energy_before=self.CalcEnergy()
            pos=(np.random.randint(0,10),np.random.randint(0,10))
            #要反转的位置 
            self.lattice[pos] = -self.lattice[pos]
            energy_after=self.CalcEnergy()
            if energy_after <=  energy_before:
                pass
                print("move")
            else:
                acc_probability= min(1.0,exp(-(energy_after - energy_before)/(k * self.temp)))
                if np.random.uniform(0, 1) <= acc_probability:
                    pass
                    print("accept")
                else: 
                    self.lattice[pos] = -self.lattice[pos]
                    print("not accept")
                    energy_after=self.CalcEnergy()#把能量算回去
            with open("move.dat","a") as f:
                f.write(str(step))
                f.write("\t")
                f.write(str(energy_after))
                f.write('\n')
            self.energy_square.append(energy_after**2)
            energy_sum += energy_after
            step+=1
            self.M_per_step.append(np.sum(self.lattice))
        #上面是每一步，下面是结束循环之后
        self.energy_average=energy_sum/self.nsteps
        self.M_average=sum(self.M_per_step)/self.nsteps

       
i2=IsingModel(10,10,0.015)
i2.MCmove()
# print(i2.M)  直接i2.M会获得0，应该i2.M=i2.Calc_M()，再打印，或者print(i2.Calc_M())
#print(i2.M_per_step)
#print(i2.M_average)
# print(i2.Calc_X())
# print(i2.X)
#后面还要再加



