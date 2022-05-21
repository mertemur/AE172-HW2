import math
import matplotlib.pyplot as plt
import numpy as np

w = 9600
s = 16
PA_max_sea = 150000
CD0 = 0.03
K0 = -0.03
K = 0.06
CL_max = 1.6
rho0 = 1.225
rhoT = [1.225,1.16727,1.11164,1.05807,1.00649,0.956859,0.909122,0.863229,0.819129,0.776775,0.736116,0.697106,0.659697,0.623844,0.589501,0.556624,0.525168,0.495090,0.466348,0.438901,0.412707,0.387725,0.363918,0.336327,0.310828,0.287262,0.265483,0.245355,0.226753,0.209562,0.193674]

def CLminforPR():
    cl = ((-math.sqrt((K0**2)+12*K*CD0))-K0) / ((-2)*K)
    return cl

def CL(v,r):
    cl = ((2*w)/(r*(v**2)*s))
    return cl

def CD(v,r):
    cd = CD0 + K0 * CL(v,r) + K * (CL(v,r)**2)
    return cd

def PowerReq(v,r):
    PR = ThrustReq(v,r) * v
    return PR

def ThrustReq(v,r):
    TR = 0.5 * r * (v**2) * s * CD(v,r)
    return TR

def PowerAv(r):
    PA = (r/rho0) * PA_max_sea
    return PA

def ExcessPower(v,r):
    EP = PowerAv(r) - PowerReq(v,r)
    return EP

def RateOfClimb(ep):
    RC = ep/w
    return RC

def PathAngleRad(v,rc):
    angle = math.asin(rc/v)
    return(angle)

def Vstall(r):
    vs = math.sqrt((2*w)/(r*s*CL_max))
    return vs

def VforminPR(r):
    vm = math.sqrt((2*w) / (r*s*CLminforPR()))
    return vm

def RhoAbsCeiling():
    rh = ((math.sqrt((2*(w**3))/s)*4*CD0*rho0)/((CL_max**(3/2))*PA_max_sea))**(2/3)
    #rh = ((CD0/(CL_max**2)+K0/CL_max+1)*(2*(w**2)*rho0)*math.sqrt((s*CL_max)/(2*w))/(s*PA_max_sea))**(2/3)
    return rh



def firstQ():
    v = VforminPR(rho0)
    print("1) ","V_EP_max = ", v, " Rate of Climb = ", RateOfClimb(ExcessPower(v,rho0)))

def secondQ():
    rho = ((PowerReq(Vstall(rho0),rho0)/PowerAv(rho0))**(2/3))*rho0
    print("2) ","From table, look for the altitude with density = ", rho)

def thirdQ():
    rac = RhoAbsCeiling()
    print("3) ","From table, look for the altitude with density = ", rac, " Speed = ", VforminPR(rac))

def fourthQ():
    rc = RateOfClimb(ExcessPower(50,rho0))
    deg = (PathAngleRad(50,rc)*180)/math.pi
    print("4) ","Rate of Climb = ", rc, " Flight Path Angle = ", deg)

def fifthQ():
    maxa = ((4*K*(w**2))/((50)*s))/(math.sqrt(((K0*w*50)**2)-8*(K*(w**2)/((50)*s))*((0.5*(50**3)*s*CD0)-(PA_max_sea/rho0)))-(K0*w*50))
    h = np.arange(0,15001,500)
    rc =[]
    for idx, i in enumerate(rhoT):
        trc = RateOfClimb(ExcessPower(50,rhoT[idx]))
        rc.append(trc)
        if (trc <= 0):
            h = np.arange(0,(500*idx+1),500)
            break
    plt.plot(h,rc)
    plt.title("Q5 - RC v h")
    plt.xlabel("h (m)")
    plt.ylabel("Rate of Climb")
    plt.grid()
    plt.savefig("Q5_RC_v_h.png", dpi = 300, bbox_inches = "tight")
    plt.show()
    print("5) ", "Q5_RC_v_h.png", " Maximum altitude (look from table) with V=50 = ", maxa)


def sixthQ():
    h = np.arange(0,15001,500)
    rc =[]
    for idx, i in enumerate(rhoT):
        trc = CL(50,i)*math.cos(PathAngleRad(50,RateOfClimb(ExcessPower(50,i))))
        rc.append(trc)
    plt.plot(h,rc)
    plt.title("Q6 - CL v h")
    plt.xlabel("h (m)")
    plt.ylabel("CL")
    plt.grid()
    plt.savefig("Q6_CL_v_h.png", dpi = 300, bbox_inches = "tight")
    plt.show()

    h = np.arange(0,15001,500)
    rc =[]
    for idx, i in enumerate(rhoT):
        trc = PathAngleRad(50,RateOfClimb(ExcessPower(50,i)))*180/math.pi
        rc.append(trc)
        if (trc <= 0):
            h = np.arange(0,(500*idx+1),500)
            break
    plt.plot(h,rc)
    plt.title("Q6 - gamma v h")
    plt.xlabel("h (m)")
    plt.ylabel("Flight Path Angle(degree)")
    plt.grid()
    plt.savefig("Q6_gamma_v_h.png", dpi = 300, bbox_inches = "tight")
    plt.show()
    print("6) ","Q6_CL_v_h.png","Q6_gamma_v_h.png")

def seventhQ():
    h = np.arange(0,15001,500)
    rc =[]
    for idx, i in enumerate(rhoT):
        trc = VforminPR(i)
        if (trc >= VforminPR(RhoAbsCeiling())):
            
            rc.append(VforminPR(RhoAbsCeiling()))
            h = np.arange(0,(500*(idx)+1),500)
            break
        rc.append(trc)
        
    plt.plot(h,rc)
    plt.title("Q7 - V_EP_max v h")
    plt.xlabel("h (m)")
    plt.ylabel("V_EP_max")
    plt.grid()
    plt.savefig("Q7V_EP_max_v_h.png", dpi = 300, bbox_inches = "tight")
    plt.show()
    print("7) ","Q7V_EP_max_v_h.png")


firstQ()
secondQ()
thirdQ()
fourthQ()
fifthQ()
sixthQ()
seventhQ()