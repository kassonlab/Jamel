import numpy as np
def read_DLS_histo(csv_file):

    data=[['Diameter','Intensity']]
    with open(csv_file,'r') as file:
        for row in file:
            values=row.split(',')
            # checking for rows with 2d info
            if len(values)==2:
                # cutting off data extraction after collecting peak intensity data with new title identifier "number'
                if values[1].startswith('Number'):
                    break
                else:
                    values=[float(x) for x in values if x[0].isdigit()]
                    if values and values[1]!=0:
                        data.append(values)
        data=np.array(data,dtype=object)
        data[1:,0]*=2
        return data
file_array=read_DLS_histo("C:\\Users\jamel\Downloads\\attachments\Bloom\\t2test-D8.noname-ma-histo.csv")
#look at range, stddev, translate to diameter
# def sum_radius_intensity(radius_data,bin_size):

print(file_array)