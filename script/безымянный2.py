import matplotlib.pyplot as plt
import csv


path_individual_1 = [
"D:/pythonProject/MitoFragility/MitoFragilityScore/Fragility/SEQ-g38_Mt-Short_Test-test_individual_1/SEQ-g38_Mt-Short_Test-test_individual_1-CGS-50-50-50-1000-500-50-CEN-50-FF.csv",
"D:/pythonProject/MitoFragility/MitoFragilityScore/Fragility/SEQ-g38_Mt-Short_Test-test_individual_1/SEQ-g38_Mt-Short_Test-test_individual_1-CGS-50-50-50-1000-500-50-CEN-100-FF.csv",
"D:/pythonProject/MitoFragility/MitoFragilityScore/Fragility/SEQ-g38_Mt-Short_Test-test_individual_1/SEQ-g38_Mt-Short_Test-test_individual_1-CGS-50-50-50-1000-500-50-CEN-150-FF.csv",
"D:/pythonProject/MitoFragility/MitoFragilityScore/Fragility/SEQ-g38_Mt-Short_Test-test_individual_1/SEQ-g38_Mt-Short_Test-test_individual_1-CGS-50-50-50-1000-500-50-CEN-200-FF.csv",
"D:/pythonProject/MitoFragility/MitoFragilityScore/Fragility/SEQ-g38_Mt-Short_Test-test_individual_1/SEQ-g38_Mt-Short_Test-test_individual_1-CGS-50-50-50-1000-500-50-CEN-250-FF.csv",
"D:/pythonProject/MitoFragility/MitoFragilityScore/Fragility/SEQ-g38_Mt-Short_Test-test_individual_1/SEQ-g38_Mt-Short_Test-test_individual_1-CGS-50-50-50-1000-500-50-CEN-300-FF.csv",
"D:/pythonProject/MitoFragility/MitoFragilityScore/Fragility/SEQ-g38_Mt-Short_Test-test_individual_1/SEQ-g38_Mt-Short_Test-test_individual_1-CGS-50-50-50-1000-500-50-CEN-350-FF.csv",
"D:/pythonProject/MitoFragility/MitoFragilityScore/Fragility/SEQ-g38_Mt-Short_Test-test_individual_1/SEQ-g38_Mt-Short_Test-test_individual_1-CGS-50-50-50-1000-500-50-CEN-400-FF.csv",
"D:/pythonProject/MitoFragility/MitoFragilityScore/Fragility/SEQ-g38_Mt-Short_Test-test_individual_1/SEQ-g38_Mt-Short_Test-test_individual_1-CGS-50-50-50-1000-500-50-CEN-450-FF.csv",
"D:/pythonProject/MitoFragility/MitoFragilityScore/Fragility/SEQ-g38_Mt-Short_Test-test_individual_1/SEQ-g38_Mt-Short_Test-test_individual_1-CGS-50-50-50-1000-500-50-CEN-500-FF.csv"
    ]
path_individual_2 = [
"D:/pythonProject/MitoFragility/MitoFragilityScore/Fragility/SEQ-g38_Mt-Short_Test-test_individual_2/SEQ-g38_Mt-Short_Test-test_individual_2-CGS-50-50-50-1000-500-50-CEN-50-FF.csv",
"D:/pythonProject/MitoFragility/MitoFragilityScore/Fragility/SEQ-g38_Mt-Short_Test-test_individual_2/SEQ-g38_Mt-Short_Test-test_individual_2-CGS-50-50-50-1000-500-50-CEN-100-FF.csv",
"D:/pythonProject/MitoFragility/MitoFragilityScore/Fragility/SEQ-g38_Mt-Short_Test-test_individual_2/SEQ-g38_Mt-Short_Test-test_individual_2-CGS-50-50-50-1000-500-50-CEN-150-FF.csv",
"D:/pythonProject/MitoFragility/MitoFragilityScore/Fragility/SEQ-g38_Mt-Short_Test-test_individual_2/SEQ-g38_Mt-Short_Test-test_individual_2-CGS-50-50-50-1000-500-50-CEN-200-FF.csv",
"D:/pythonProject/MitoFragility/MitoFragilityScore/Fragility/SEQ-g38_Mt-Short_Test-test_individual_2/SEQ-g38_Mt-Short_Test-test_individual_2-CGS-50-50-50-1000-500-50-CEN-250-FF.csv",
"D:/pythonProject/MitoFragility/MitoFragilityScore/Fragility/SEQ-g38_Mt-Short_Test-test_individual_2/SEQ-g38_Mt-Short_Test-test_individual_2-CGS-50-50-50-1000-500-50-CEN-300-FF.csv",
"D:/pythonProject/MitoFragility/MitoFragilityScore/Fragility/SEQ-g38_Mt-Short_Test-test_individual_2/SEQ-g38_Mt-Short_Test-test_individual_2-CGS-50-50-50-1000-500-50-CEN-350-FF.csv",
"D:/pythonProject/MitoFragility/MitoFragilityScore/Fragility/SEQ-g38_Mt-Short_Test-test_individual_2/SEQ-g38_Mt-Short_Test-test_individual_2-CGS-50-50-50-1000-500-50-CEN-400-FF.csv",
"D:/pythonProject/MitoFragility/MitoFragilityScore/Fragility/SEQ-g38_Mt-Short_Test-test_individual_2/SEQ-g38_Mt-Short_Test-test_individual_2-CGS-50-50-50-1000-500-50-CEN-450-FF.csv",
"D:/pythonProject/MitoFragility/MitoFragilityScore/Fragility/SEQ-g38_Mt-Short_Test-test_individual_2/SEQ-g38_Mt-Short_Test-test_individual_2-CGS-50-50-50-1000-500-50-CEN-500-FF.csv"
    ]
path_individual_3 = [
"D:/pythonProject/MitoFragility/MitoFragilityScore/Fragility/SEQ-g38_Mt-Short_Test-test_individual_3/SEQ-g38_Mt-Short_Test-test_individual_3-CGS-50-50-50-1000-500-50-CEN-50-FF.csv",
"D:/pythonProject/MitoFragility/MitoFragilityScore/Fragility/SEQ-g38_Mt-Short_Test-test_individual_3/SEQ-g38_Mt-Short_Test-test_individual_3-CGS-50-50-50-1000-500-50-CEN-100-FF.csv",
"D:/pythonProject/MitoFragility/MitoFragilityScore/Fragility/SEQ-g38_Mt-Short_Test-test_individual_3/SEQ-g38_Mt-Short_Test-test_individual_3-CGS-50-50-50-1000-500-50-CEN-150-FF.csv",
"D:/pythonProject/MitoFragility/MitoFragilityScore/Fragility/SEQ-g38_Mt-Short_Test-test_individual_3/SEQ-g38_Mt-Short_Test-test_individual_3-CGS-50-50-50-1000-500-50-CEN-200-FF.csv",
"D:/pythonProject/MitoFragility/MitoFragilityScore/Fragility/SEQ-g38_Mt-Short_Test-test_individual_3/SEQ-g38_Mt-Short_Test-test_individual_3-CGS-50-50-50-1000-500-50-CEN-250-FF.csv",
"D:/pythonProject/MitoFragility/MitoFragilityScore/Fragility/SEQ-g38_Mt-Short_Test-test_individual_3/SEQ-g38_Mt-Short_Test-test_individual_3-CGS-50-50-50-1000-500-50-CEN-300-FF.csv",
"D:/pythonProject/MitoFragility/MitoFragilityScore/Fragility/SEQ-g38_Mt-Short_Test-test_individual_3/SEQ-g38_Mt-Short_Test-test_individual_3-CGS-50-50-50-1000-500-50-CEN-350-FF.csv",
"D:/pythonProject/MitoFragility/MitoFragilityScore/Fragility/SEQ-g38_Mt-Short_Test-test_individual_3/SEQ-g38_Mt-Short_Test-test_individual_3-CGS-50-50-50-1000-500-50-CEN-400-FF.csv",
"D:/pythonProject/MitoFragility/MitoFragilityScore/Fragility/SEQ-g38_Mt-Short_Test-test_individual_3/SEQ-g38_Mt-Short_Test-test_individual_3-CGS-50-50-50-1000-500-50-CEN-450-FF.csv",
"D:/pythonProject/MitoFragility/MitoFragilityScore/Fragility/SEQ-g38_Mt-Short_Test-test_individual_3/SEQ-g38_Mt-Short_Test-test_individual_3-CGS-50-50-50-1000-500-50-CEN-500-FF.csv",
    ]

individuals = [path_individual_1, path_individual_2, path_individual_3]
X = []
Y = []

for i in individuals:
    
    for p in i:
        X = []
        Y = []
        with open(p, 'r') as datafile:
            plotting = csv.reader(datafile, delimiter=',')
            first_row = next(plotting)
    
            for ROWS in plotting:
    
    
                Y.append(float(ROWS[-2]))
        
        plt.plot(Y)
        plt.title(f'{p.split("/")[-1]}')
        plt.xlabel('Window')
        plt.ylabel(f'{first_row[-2]}')
        plt.show()
