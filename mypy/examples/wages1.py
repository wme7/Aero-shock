'''if-else example:  calculate weekly wages -- alternate version'''

def calcWeeklyWages(totalHours, hourlyWage):  #NEW
    '''Return the total weekly wages for a worker working totalHours,
    with a given regular hourlyWage.  Include overtime for hours over 40.
    '''
    if totalHours <= 40:
        regularHours = totalHours
        overtime = 0
    else:
        overtime = totalHours - 40
        regularHours = 40
    return hourlyWage*regularHours + (1.5*hourlyWage)*overtime

def main():
    hours = float(input('Enter hours worked: '))
    wage = float(input('Enter dollars paid per hour: '))
    total = calcWeeklyWages(hours, wage)
    print('Wages for {hours} hours at ${wage:.2f} per hour are ${total:.2f}.'
          .format(**locals()))

main()
