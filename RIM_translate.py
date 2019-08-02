def RIM_translate(rim:str):
    '''This function takes in the RIM matrix, and if applicable returns the
    known RIM code as defined in the RIM paper; otherwise returns None'''
    rims_in_paper = {
        '[0.5, 0.5, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0]' : 'RIM 1',
        '[0.5, 0.5, 1.0, 0.0, 0.5, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0]' : 'RIM 2',
        '[0.5, 0.5, 0.5, 0.0, 0.5, 0.5, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0]' : 'RIM 3',
        '[0.5, 0.5, 1.0, 0.0, 0.0, 1.0, 0.0, 0.5, 1.0, 0.0, 0.0, 1.0]' : 'RIM 4',
        '[0.5, 0.5, 1.0, 0.0, 0.5, 1.0, 0.0, 0.5, 1.0, 0.0, 0.5, 1.0]' : 'RIM 5',
        '[0.5, 0.5, 0.5, 0.0, 0.5, 0.5, 0.0, 0.5, 1.0, 0.0, 0.5, 1.0]' : 'RIM 6',
        '[0.5, 0.5, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0]' : 'RIM 7',
        '[0.5, 0.5, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0]' : 'RIM 8',
        '[0.5, 0.5, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0]' : 'RIM 9',
        '[0.5, 0.5, 1.0, 0.0, 0.0, 1.0, 0.5, 1.0, 1.0, 0.0, 0.0, 1.0]' : 'RIM 10',
        '[0.5, 1.0, 1.0, 0.0, 1.0, 1.0, 0.5, 1.0, 1.0, 0.0, 1.0, 1.0]' : 'RIM 11',
        '[0.5, 0.5, 0.0, 0.0, 1.0, 0.0, 0.5, 0.5, 1.0, 0.0, 1.0, 1.0]' : 'RIM 12',
        '[1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0]' : 'RIM 13',
        '[1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0]' : 'RIM 14',
        '[1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0]' : 'RIM 15',
        '[0.5, 0.5, 1.0, 0.0, 0.0, 1.0, 0.5, 0.5, 1.0, 0.0, 0.0, 1.0]' : 'RIM 16',
        '[0.5, 0.5, 1.0, 0.0, 0.5, 1.0, 0.5, 0.5, 1.0, 0.0, 0.5, 1.0]' : 'RIM 17',
        '[0.5, 0.5, 0.5, 0.0, 0.5, 0.5, 0.5, 0.0, 0.5, 0.0, 0.5, 0.5]' : 'RIM 18',
        '[0.0, 0.5, 1.0, 0.0, 0.0, 1.0, 0.0, 0.5, 1.0, 0.0, 0.0, 1.0]' : 'RIM 19',
        '[0.0, 0.5, 1.0, 0.0, 0.5, 1.0, 0.0, 0.5, 1.0, 0.0, 0.5, 1.0]' : 'RIM 20',
        '[0.0, 0.5, 0.5, 0.0, 0.5, 0.5, 0.0, 0.5, 0.5, 0.0, 0.5, 0.5]' : 'RIM 21',
        '[0.0, 0.0, 1.0, 0.0, 0.5, 1.0, 0.0, 0.0, 1.0, 0.0, 0.5, 1.0]' : 'RIM 22',
        '[0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0]' : 'RIM 23',
        '[0.5, 0.5, 1.0, 0.0, 0.5, 1.0, 0.5, 0.5, 1.0, 0.0, 0.0, 1.0]' : 'RIM 24',
        '[1.0, 1.0, 1.0, 0.0, 0.5, 1.0, 1.0, 1.0, 1.0, 0.0, 0.5, 1.0]' : 'RIM 25',
        '[1.0, 1.0, 1.0, 0.0, 0.5, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0]' : 'RIM 26',
        '[0.5, 0.5, 0.5, 0.0, 0.5, 0.5, 0.5, 0.5, 1.0, 0.0, 0.5, 1.0]' : 'RIM 27',
        '[0.5, 0.5, 0.5, 0.0, 0.5, 0.5, 0.5, 0.5, 1.0, 0.0, 0.0, 1.0]' : 'RIM 28',
    }
    translated_rim = rims_in_paper.get(rim)
    return translated_rim