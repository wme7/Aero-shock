'''Compare print with concatenation and with two versions of format string.'''

applicant = input("Enter the applicant's name: ")
interviewer = input("Enter the interviewer's name: ")
time = input("Enter the appointment time: ")
print(interviewer + ' will interview ' + applicant + ' at ' + time +'.')
print(interviewer, ' will interview ', applicant, ' at ', time, '.', sep='')
print('{} will interview {} at {}.'.format(interviewer, applicant, time))
print('{0} will interview {1} at {2}.'.format(interviewer, applicant, time))
