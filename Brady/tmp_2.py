
def numericable_str(my_string):

    is_numericable_str = False
    try:
        my_numeric_value = int(my_string)
        is_numericable_str = True
    except ValueError:
        try:
            my_numeric_value = float(my_string)
            is_numericable_str = True
        except ValueError:
            is_numericable_str = False

    return is_numericable_str


