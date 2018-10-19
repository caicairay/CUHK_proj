import sys

def dictionary(amount_in, unit_in, unit_out):
    length_dict = {'cm':1, 'm':0.01, 'km':1e-5, 'pc':3.2407793e-19}
    mass_dict = {'g':1, 'kg':0.001, 'msun':5.0287e-34}
    time_dict = {'s':1, 'myr':3.168876e-14}
    density_dict = {'g/cm3':1, 'msun/pc3':1.477432e+22}
    bfield_dict = {'muG':1,'G':1e-6,'code_unit':1.2430766}
    
    # also check the unit_out, if they are not in the same dict, return an error message
    if unit_in in length_dict.keys():
      return float(amount_in)*length_dict[unit_out]/length_dict[unit_in]
    
    elif unit_in in mass_dict.keys():
      return float(amount_in)*mass_dict[unit_out]/mass_dict[unit_in]

    elif unit_in in time_dict.keys():
      return float(amount_in)*time_dict[unit_out]/time_dict[unit_in]  

    elif unit_in in density_dict.keys():
      return float(amount_in)*density_dict[unit_out]/density_dict[unit_in]

    elif unit_in in bfield_dict.keys():
      return float(amount_in)*bfield_dict[unit_out]/bfield_dict[unit_in] 
      
    else: return 0   
#The main function of this unit converter
def main():
  amount_in, unit_in, unit_out = sys.argv[1], sys.argv[2], sys.argv[3]
  amount_out = dictionary(amount_in, unit_in, unit_out)
  if bool(amount_out): 
    #print amount_in, unit_in, '=', amount_out, unit_out
    print amount_out
    #print round(amount_out,6)
  else:
    print 'Unit not found!'
if __name__ == "__main__":
  main()
