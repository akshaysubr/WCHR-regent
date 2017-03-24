import matplotlib.pyplot as plt
import numpy
import h5py
import sys

def get_regent_TKE(data_file):
    t   = []
    TKE = []
    
    f = open(data_file)
    
    for line in f:
        line = line.strip().split()
        if line and line[0].isdigit():
            t.append(float(line[2]))
            TKE.append(float(line[12]))
    
    f.close()
    
    t = numpy.array(t)
    TKE = numpy.array(TKE)

    return t, TKE

# if len(sys.argv) != 3:
#   print "Usage:"
#   print "  python %s <data file> <linetype>" % sys.argv[0]
#   sys.exit()
# 
# data_file  = sys.argv[1]
# linetype   = sys.argv[2]

t_ECISFD, TKE_ECISFD = get_regent_TKE('../src/taylorgreen_ECISFD_p1_N032.log')
t_ECIEFD, TKE_ECIEFD = get_regent_TKE('../src/taylorgreen_ECIEFD_p1_N032.log')
t_ECIEFD2, TKE_ECIEFD2 = get_regent_TKE('../src/taylorgreen_ECIEFD_N032.log')
t_CI, TKE_CI = get_regent_TKE('../src/taylorgreen_N032.log')
t_EI, TKE_EI = get_regent_TKE('../src/taylorgreen_LD_N032.log')

t_JS, TKE_JS = numpy.loadtxt('taylorgreen_results/32x32x32/WCNS5-JS_KE.txt', unpack=True)
t_Z,  TKE_Z  = numpy.loadtxt('taylorgreen_results/32x32x32/WCNS5-Z_KE.txt', unpack=True)
t_CU, TKE_CU = numpy.loadtxt('taylorgreen_results/32x32x32/WCNS6-CU-M2_KE.txt', unpack=True)
t_LD, TKE_LD = numpy.loadtxt('taylorgreen_results/32x32x32/WCNS6-LD_KE.txt', unpack=True)

plt.figure(1)
plt.plot(t_CI,     TKE_CI,     'b-', linewidth=1,fillstyle='none',label=r'$WCHR6-CI$')
plt.plot(t_ECIEFD2, TKE_ECIEFD2, 'k-.', linewidth=2,fillstyle='none',label=r'$WCHR6-ECI-EFD \; p=2$')
plt.plot(t_ECIEFD, TKE_ECIEFD, 'k--', linewidth=2,fillstyle='none',label=r'$WCHR6-ECI-EFD \; p=1$')
plt.plot(t_ECISFD, TKE_ECISFD, 'k-', linewidth=1,fillstyle='none',label=r'$WCHR6-ECI-SFD \; p=1$')
# plt.plot(t_JS,TKE_JS/TKE_JS[0],'c--',linewidth=1,fillstyle='none',label=r'$WCNS5-JS$')
plt.plot(t_Z, TKE_Z /TKE_Z [0],'r-.',linewidth=1,fillstyle='none',label=r'$WCNS5-Z $')
# plt.plot(t_CU,TKE_CU/TKE_CU[0],'m-.',linewidth=2,fillstyle='none',label=r'$WCNS6-CU-M2$')
plt.plot(t_LD,TKE_LD/TKE_LD[0],'g--',linewidth=2,fillstyle='none',label=r'$WCNS6-LD$')
plt.plot(t_EI,TKE_EI,          'm-.',linewidth=2,fillstyle='none',label=r'$WCNS6-EI$')

plt.xlim((0,10))
plt.xlabel(r'$t$', fontsize=20)
plt.ylabel(r'$TKE$', fontsize=20)
plt.legend(loc='lower left')
plt.show()
