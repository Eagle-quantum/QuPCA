import warnings

class customWarning():
    
    """
    Class to customize the warning messages. It override the warn method reporting only the desired warning message.

        Parameters
        ----------
        
        Notes
        ----------
    """
    
    
    def custom_formatwarning(msg, *args, **kwargs):
        # ignore everything except the message
        return str(msg) + '\n'
    
    warnings.formatwarning = custom_formatwarning
    
    @classmethod
    def warn(cls,message):
        warnings.warn(message)