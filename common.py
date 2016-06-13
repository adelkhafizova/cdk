smf_active = set()
smf_inactive = set()
smf_inactive_wrong = set()
dragon_active = set()
dragon_inactive = set()
dragon_inactive_wrong = set()
full_set_active = set()
full_set_inactive = set()
full_set_inactive_wrong = set()

with open('cyp_smf_statistics.csv', 'r') as f:
    smf_buffer = f.read()
with open('cyp_dragon_statistics.csv', 'r') as f:
    dragon_buffer = f.read()
with open('cyp_full_set_statistics.csv', 'r') as f:
    full_set_buffer = f.read()
    
smf_data = list()
dragon_data = list()
full_set_data = list()

smf_buffer = smf_buffer.split('\n')
for line in smf_buffer:
    smf_data.append(line.split('\t'))

dragon_buffer = dragon_buffer.split('\n')
for line in dragon_buffer:
    dragon_data.append(line.split('\t'))

full_set_buffer = full_set_buffer.split('\n')
for line in full_set_buffer:
    full_set_data.append(line.split('\t'))

for properties in smf_data:
    try:
        if float(properties[-1]) > 0.85 and int(properties[-2]) == 6:
            smf_active.add(properties[0])
    except:
        pass

for properties in dragon_data:
    try:
        if float(properties[-1]) > 0.85 and int(properties[-2]) == 6:
            dragon_active.add(properties[0])
    except:
        pass

for properties in full_set_data:
    try:
        if float(properties[-1]) > 0.85 and int(properties[-2]) == 6:
            full_set_active.add(properties[0])
    except:
        pass

for properties in smf_data:
    try:
        if float(properties[-1]) > 0.85 and int(properties[-2]) == 7:
            smf_inactive.add(properties[0])
    except:
        pass

for properties in dragon_data:
    try:
        if float(properties[-1]) > 0.85 and int(properties[-2]) == 7:
            dragon_inactive.add(properties[0])
    except:
        pass

for properties in full_set_data:
    try:
        if float(properties[-1]) > 0.85 and int(properties[-2]) == 7:
            full_set_inactive.add(properties[0])
    except:
        pass
        
for properties in smf_data:
    try:
        if float(properties[-1]) > 0.85 and int(properties[-2]) == 4:
            smf_inactive_wrong.add(properties[0])
    except:
        pass

for properties in dragon_data:
    try:
        if float(properties[-1]) > 0.85 and int(properties[-2]) == 4:
            dragon_inactive_wrong.add(properties[0])
    except:
        pass

for properties in full_set_data:
    try:
        if float(properties[-1]) > 0.85 and int(properties[-2]) == 4:
            full_set_inactive_wrong.add(properties[0])
    except:
        pass
        
active = smf_active & dragon_active & full_set_active
inactive = smf_inactive & dragon_inactive & full_set_inactive
active_to_inactive = smf_inactive_wrong & dragon_inactive_wrong & full_set_inactive_wrong

print active
print inactive
print active_to_inactive
