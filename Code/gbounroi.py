# gbounroi python
# the following is a python version of the matlab gbounroi
# written in matlab by TBB 2012
# Translation to py written by AHB 2020
# These functions use existing PyQt5 libraries to create a GUI


def gbounroi(*varargin):
    gui_singleton = 1
    gui_state = dict(
        gui_Name = filename
        guisingleton = gui_singleton
        gui_openingfcn = @gbounroi_openingfcn
        gui_outputfcn = @gbounroi_outputfcn
        gui_layoutfcn = []
        gui_callback = []
    )

    if varagin and str(varargin[0]):
        # should convert first arg to function handle
        gui_state[gui_callback] = varargin[0] 


    if nargout:
        # if argout is true go to guimainfunct
        pass
    else:
        pass
    
    return varargout

def gbounroi_openingfcn(h):
    pass

def gbounroi_outputfcn():
    pass

