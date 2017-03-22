# The Insertion Operation in Error Generation

The whole bunch of stuff is implemented in a long loop.

#### 1. Initialization for creating duplicate record(s). 

Febrl does not rely on data in real world, so this step is equivalent with interation through all the records.

`max_num_record_modifi` is set outside the loop, which restricts the maximum number of errors allowed in one record.

```python
num_modif_in_record = 0  # Count the number of modifications in this record
field_mod_count_dict = {} # Set the field modification counters to zero for all fields
```
#### 2. While `num_modif_in_record` smaller than `max_num_record_modifi`, do modification
```python
# Randomly choose a field
field_dict = random_select(select_prob_list)
field_name = field_dict['name']

# Make sure this field hasn't been modified already
while (field_mod_count_dict[field_name] == max_num_field_modifi):
    field_dict = random_select(select_prob_list)
    field_name = field_dict['name']
```
Define the valid values for each field (This is important for the insertion operation):
```python
if (field_dict['char_range'] == 'digit'):
    field_range = string.digits
elif (field_dict['char_range'] == 'alpha'):
    field_range = string.lowercase
elif (field_dict['char_range'] == 'alphanum'):
    field_range = string.digits+string.lowercase
```
Make sure the number of modifications do not exceed the maximum restriction:
```python
if (max_num_field_modifi == 1):
    num_field_mod_to_do = 1
else:
    num_field_mod_to_do = random.randint(1, max_num_field_modifi)
num_rec_mod_to_do = max_num_record_modifi - num_modif_in_record
```
Initialization for the modifications in the field:
```python
num_modif_in_field = 0  # Count the number of modifications in this field

org_field_val = org_rec_dict.get(field_name, None) # Get original value
```
#### 3. Loop over chosen number of modifications
Randomly choose a modification. In Febrl, each field has different probabilities to be chosen, and so does each modification.
The probabilities are set outside the loop.`field_dict['prob_list']` provides all the possible modification 
as well as their probabilities.
```python
mod_op = random_select(field_dict['prob_list'])
```
Then do the modification. If it's an insertion modification and the value in the field is not empty,
choose position and character of insertion.
```python
rand_ins_pos = error_position(dup_field_val, +1)
rand_char = random.choice(field_range)
```
Insertion:
```python
if (rand_ins_pos != None):  # If a valid position was returned
  dup_field_val = dup_field_val[:rand_ins_pos] + rand_char + \
                  dup_field_val[rand_ins_pos:]

  if (VERBOSE_OUTPUT == True):
    print '      Inserted char "%s" into field "%s": "%s" -> "%s"' \
          % (rand_char, field_name, old_field_val, dup_field_val)
```
#### 4. Check if the modification is successful, and update the modification counts.
A modification is successful if it makes a difference to the record than the original record, and it is not a duplicate with 
another record.

Update the field modification counter and record modification counter:
```python
if (old_field_val == org_field_val) and \
   (dup_field_val != old_field_val):  # The first field modification
  field_mod_count_dict[field_name] = 1
  num_modif_in_record += 1

elif (old_field_val != org_field_val) and \
     (dup_field_val != old_field_val):  # Following field modifications
  field_mod_count_dict[field_name] += 1
  num_modif_in_record += 1

if (dup_field_val != old_field_val):
  dup_rec_dict[field_name] = dup_field_val
```
Finally, update the data set:
```python
rec_data = dup_rec_dict.copy()  # Make a copy of the record dictionary
del(rec_data['rec_id'])  # Remove the record identifier
rec_list = rec_data.items()
rec_list.sort()
rec_str = str(rec_list)

if (rec_str not in all_rec_set):  # Check if same record already created
  all_rec_set.add(rec_str)
  org_rec_used[org_rec_id] = 1

  dup_rec[dup_rec_id] = dup_rec_dict  # Insert into duplicate records
  rec_cnt += 1

  d += 1  # Duplicate counter (loop counter)
```
