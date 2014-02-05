"""
String Substitution for a Mad Lib
Adapted from code by Kirby Urner
"""                                                  

storyFormat = """                                       
Once upon a time, deep in an ancient jungle,
there lived a {animal}.  This {animal}
liked to eat {food}, but the jungle had
very little {food} to offer.  One day, an
explorer found the {animal} and discovered
it liked {food}.  The explorer took the
{animal} back to {city}, where it could
eat as much {food} as it wanted.  However,
the {animal} became homesick, so the
explorer brought it back to the jungle,
leaving a large supply of {food}.

The End
"""                                                 

def tellStory():                                     
    userPicks = dict()                              
    addPick('animal', userPicks)            
    addPick('food', userPicks)            
    addPick('city', userPicks)            
    story = storyFormat.format(**userPicks)
    print(story)
                                                    
def addPick(cue, dictionary):
    '''Prompt for a user response using the cue string,
    and place the cue-response pair in the dictionary.
    '''
    prompt = 'Enter an example for ' + cue + ': '
    response = input(prompt).strip() # 3.2 Windows bug fix
    dictionary[cue] = response                                                             

tellStory()                                         
input("Press Enter to end the program.")        
