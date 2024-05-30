import numpy as np
import os

# ORIGIN -0.5 -100.5 0
# SPACING 0.5 0.5 0.5
#ORIGIN
yo = -100.5
xo = -0.5

dx = 0.5
dy = 0.5
dim_x = 303 
dim_y = 503

xmin = 27.5
xmax = 100.0 
ymin = -0.5
ymax =  0.5

y_sal = 100 #mm

part_count = 89929

def find_integer_substring_indices(string):
    start_index = None
    end_index = None

    for i in range(len(string)):
        if string[i].isdigit():  # If the character is a digit
            if start_index is None:  # If this is the first digit encountered
                start_index = i
            end_index = i  # Update end_index on each iteration as long as digit characters are consecutive
        elif start_index is not None:  # If the consecutive digit sequence ends
            break

    return start_index

def find_integer_substring_end_indices(string):
    start_index = end_index = None

    for i in range(len(string)):
        j = len(string)-1-i
        if string[j].isdigit():  # If the character is a digit
            if end_index is None:  # If this is the first digit encountered
                end_index = j
            # print ("END DIG " + string[j])
        elif end_index is not None:  # If the consecutive digit sequence ends
            break
    
    #NOW first number
    end = False
    j = end_index - 1
    while (not end):
      if (not string[j].isdigit()):
        start_index  = j + 1
        end = True
      j -=1
    return start_index, end_index

    
def convert_substring_to_integer(string, start_index, end_index):
    # Extract the substring
    # print("STRING %s" %string)
    substring = string[start_index:end_index+1]
    print("SUBS -%s-" %substring)
    # # Convert the substring to an integer
    try:
        integer_value = int(substring)
        print ("int %d" %integer_value)
        return integer_value
    except ValueError:
        print("Error: The substring cannot be converted to an integer.")
        return None


def write_list(force):
  fi_x = open("force.csv","w")
  fi_x.write("t,f\n")

  dt = 1.0e-3
  t = 0.0
  for i in range(len(force)):
    fi_x.write(str(t) + ", " + str(force [i]) + "\n" )
    print(force [i])
    t +=dt

def open_files_with_extension(directory, extension):
    # Check if the directory exists
    if not os.path.isdir(directory):
        print(f"Directory '{directory}' does not exist.")
        return
    
    # Get a list of files with the specified extension
    files = [f for f in os.listdir(directory) if f.endswith(extension)]
    
    if not files:
        print(f"No files found with '{extension}' extension in '{directory}'.")
        return
    
    
    print("Found %d files " %len(files))
    force = []
    force = np.zeros(len(files))


    #GENERATING LIST------
    idx_list = []
    for file_name in files:
      start_idx, end_idx = find_integer_substring_end_indices(file_name)
      print("start: " + str(start_idx) + ", end " + str(end_idx))
      file_path = os.path.join(directory, file_name)

      with open(file_path, 'r') as f:
        print("File " + file_name + " found")
        idx = convert_substring_to_integer(file_name,start_idx,end_idx)
        print ("Index ", idx)
        idx_list.append(idx)
    
    print ("index list ", idx_list)
    idx_sort = sorted(range(len(idx_list)), key=lambda k: idx_list[k])
    
    print ("sorted idx", idx_sort)
    
    first = True 
    # Open each file
    tot = len(files)
    file_i = 0
    for file_name in files:
      start_idx, end_idx = find_integer_substring_end_indices(file_name)
      print("start: " + str(start_idx) + ", end " + str(end_idx))
      file_path = os.path.join(directory, file_name)


  # try:
      with open(file_path, 'r') as f:
        print("File " + file_name + " found")
        # idx = convert_substring_to_integer(file_name,start_idx,end_idx)
        idx = idx_sort[file_i]
        # print ("idx %d" %idx) 
        lines = f.readlines()
        print ("file " + str(file_i) + "/"+ str(tot))
        if (first):  
          i = 0 
          end = False  
          while (not end):              
            if (lines[i].find("VECTORS contactforce double") != -1):
              print("FOUND STRING. in line ", i, lines[i])
              line_cf = i + 1 
              end = True    
            i = i+1
              
          print ("Found in line ", i)
          first = False
		
        
        jocount = int((ymin - yo)/dy) #count to initial point 
        iocount = int((xmin - xo)/dx)
        icount  = int((int)(xmax - xmin)/dx)
        print ("icount ", icount)
        l_o = jocount * dim_x #line to initial point
        
        i_ini = i

        y = ymin -dy
        line = line_cf + l_o #INITIAL TO Y DIE COORDINATE
        jinc = 0
        print ("Reading force ")
        tot_f = 0.0
        #while (y < ymax + dy):
        while (i < i_ini + part_count-1):
          print ("y = ", y, "----------------------")
          fy = 0.0
          for i in range(icount):
            # print ("xi ", i )
            # print ("index ", line + jinc * dim_x + iocount + i)
            # print ("LINES ", lines[line + jinc * dim_x + iocount + i])
            
            numbers_array = [float(num) for num in lines[line + jinc * dim_x + iocount + i].split()]
            # print ("convlist ",numbers_array[1])
            fy +=  numbers_array[1]
            i+=1;
          
          print ("Integrated Force ", fy)
          print ("reached line ", line + jinc * dim_x + iocount + i)
          tot_f += fy
          y += dy
          jinc += 1
          
          
        # sprint ("Force", numbers_array[0], numbers_array[1],numbers_array[2])
        #force.append(numbers_array[2])
        force[idx] =  tot_f
      #wth file
      file_i += 1
      
  # except Exception as e:
      # print(f"Error opening {file_name}: {e}")
    return force
          
###### INPUT PARAM ENTRADA ##############

# directory = '/path/to/directory'
directory = '.'
extension = '.vtk'  # Change this to the extension you want to search for

force = open_files_with_extension(directory, extension)
write_list(force)
