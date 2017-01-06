#!/home/mwebb/anaconda/bin/python

#==================================================================
#  IMPORT MODULES
#==================================================================
from math import *
import sys,argparse
from numpy import *
import glob
from scipy.optimize import curve_fit

#==================================================================
#  AUX: create_parser
#==================================================================
def create_parser():

  parser = argparse.ArgumentParser(description='Computes time constants and Kohlrausch-William-Watts (kww) fits with stretching coefficients from autocorrelation decays.')
  

  # DO REQUIRED POSITIONAL ARGUMENTS
  parser.add_argument('data_name',help = 'Name of file containing data. First column should be time variable, and \
        subsequent columns should contain the ACF data for given Rouse modes.')

  # OPTIONAL ARGUMENTS
  parser.add_argument('-tcut'     ,dest='tcut',default=1.0e8,
                      help = 'Time cutoff to use for fitting data. This can be set to avoid fitting \
                              difficult long-time tails, for example. (default = 1e8)')

  parser.add_argument('-maxtau'     ,dest='tau_max',default=1.0e8,
                      help = 'Specification for constraint on optimization of time constants. \
                              If time constants beyond this value are encountered, a hyperbolic penalty is imposed. (default = 1e8)')

  parser.add_argument('-penalty'     ,dest='pen_params',default="1000 1",
                      help = 'Specification of parameters for imposing hyperbolic penalty function of the form: \
                              C{1+exp[-2k(x-xmax)]}^(-1) \
                              Should be supplied as a quoted pair value as "C k". (default = "1000 1") ')

  return parser

#==================================================================
#  AUX: convert_args
#==================================================================
def convert_args(args):
  data_name       = args.data_name
  tcut            = float(args.tcut)
  maxtau          = float(args.tau_max)
  penalty         = [float(i) for i in args.pen_params.split()]
  return(data_name,tcut,maxtau,penalty)

#==================================================================
#  AUX: is_number
#==================================================================
def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

#==================================================================
#  AUX: process_data
#==================================================================
def process_data(file):
  fid = open(file,"r")
  lines = [line.strip() for line in fid]
  lines = [line.split() for line in lines]
  lines = [line for line in lines if is_number(line[0]) ]  # SKIP HEADER LINES
  fid.close()
  x     = [float(line[0]) for line in lines if float(line[0]) > 0.0]
  y     = [[float(i) for i in line[1:]] for line in lines if float(line[0]) > 0.0]

  return (array(x),array(y).transpose())

#==================================================================
#  AUX: rouse_fun
#==================================================================
def rouse_fun(t,tau,tau_max,Cpen):
   penalty = 0.0
   if (tau >= tau_max):
     penalty += Cpen[0]/(1 + exp(-2*Cpen[1]*(tau-tau_max)))
   return exp(-t/tau) + penalty

#==================================================================
#  AUX: kww_fun
#==================================================================
def kww_fun(t,tau,beta,tau_max,Cpen):
   penalty = 0.0
   if (tau >= tau_max):
     penalty += Cpen[0]/(1 + exp(-2*Cpen[1]*(tau-tau_max)))
   if (beta < 0 ):
     penalty += Cpen[0]/(1 + exp(-2*Cpen[1]*(-beta)))
     return penalty
   return exp(-(t/tau)**beta) + penalty

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#  MAIN: _main
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
def main(argv):
  # CREATE THE ARGUMENT PARSER
  parser = create_parser()

  # PARSE ARGUMENTS
  args               = parser.parse_args()
  (data_name,tcut,taumax,Cpen) \
                     = convert_args(args)

  # PROCESS FILE 
  (t,ACF) = process_data(data_name)   # extracts data in aggregate
  N = ACF.shape[0]                    # number of data sets to be fit
  ind = t <= tcut                     # yields index reference that satisfy time cutoff
  
  # CREATE ANONYMOUS FUNCTIONS THAT IMPOSE PENALTIES
  rouse_con_fun = lambda t,tau     : rouse_fun(t,tau,taumax,Cpen)
  kww_con_fun   = lambda t,tau,beta: kww_fun(t,tau,beta,taumax,Cpen)
  
  # COMPUTE RELAXATION TIMES AND STRETCHING COEFFICIENTS
  # Fit to a functional form of the kind KWW = exp[-(t/tau)^beta]
  # and to a the rouse expected from of  R   = exp[-(t/tau)]
  tau_rouse = [0.0]*N
  tau_kww   = [0.0]*N
  beta_kww  = [0.0]*N
  ACF_kww   = zeros(ACF.shape)
  ACF_rouse = zeros(ACF.shape)
  for p,ACFp in enumerate(ACF):
    # PERFORM THE OPTIMIZATION AND CONSTRUCT FITS
    # ROUSE:
    params_rouse   = curve_fit(rouse_con_fun,t[ind],ACFp[ind],p0=[1.0],maxfev=5000)
    [tau_rouse[p]] = params_rouse[0]
    ACF_rouse[p,:] = rouse_fun(t,tau_rouse[p],1e8,[0,0])
    # kww:
    params_kww               = curve_fit(kww_con_fun,t[ind],ACFp[ind],p0=[1.0,1],maxfev=50000)
    [tau_kww[p],beta_kww[p]] = params_kww[0]
    ACF_kww[p,:]             = kww_fun(t,tau_kww[p],beta_kww[p],1e8,[0,0])
    

  # WRITE OUT FITS
  # ROUSE:
  with open("fit.exp." +data_name,"w") as fid:
    fid.write("#{:<6s}".format("t"))
    for p in range(N):
      fid.write("p={:<8d}".format(p))
    fid.write("\n")
    for i,ti in enumerate(t):
      fid.write("{:>7.3f}".format(ti))
      for p in range(N):
        fid.write("{:>10.5f}".format(ACF_rouse[p,i]))
      fid.write("\n")
  # kww:
  with open("fit.kww." +data_name,"w") as fid:
    fid.write("#{:<6s}".format("t"))
    for p in range(N):
      fid.write("p={:<8d}".format(p))
    fid.write("\n")
    for i,ti in enumerate(t):
      fid.write("{:>7.3f}".format(ti))
      for p in range(N):
        fid.write("{:>10.5f}".format(ACF_kww[p,i]))
      fid.write("\n")

  with open("taus." + data_name,"w") as fid:
    fid.write("{:<15s}{:^15s}{:^15s}{:^15s}\n".format("#p","tau_rouse","tau_kww","beta_kww"))
    for p,(tr,tkww,bff) in enumerate(zip(tau_rouse,tau_kww,beta_kww)):
      fid.write("{:<15d}{:>15.5f}{:>15.5f}{:>15.5f}\n".format(p,tr,tkww,bff))
  
  
#==================================================================
#  RUN PROGRAM
#==================================================================
if __name__ == "__main__":
  main(sys.argv[1:])
                          
