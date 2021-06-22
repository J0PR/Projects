import pythoncom, pyHook
import serial
import datetime
from time import sleep

class Keystroke_Watcher(object):
    def __init__(self):
        self.hm = pyHook.HookManager()
        self.hm.KeyDown = self.OnKeyboardEvent
        self.hm.HookKeyboard()
        self.serevent = serial.Serial()
        self.serevent.baudrate = 9600
        self.serevent.port = 'COM12'
        self.serevent.close() #make sure its closed
        self.serevent.open()
        self.filename = "V1_LOG_"+datetime.datetime.now().strftime('%Y-%m-%d H%Hm%mS%S') + ".txt"
        self.file= open(self.filename,"w+")
        self.file.close()


    def OnKeyboardEvent(self, event):
        self.serevent.write(event.KeyID.to_bytes(1,'little'))
        sleep(0.01)
        aux=0
        self.serevent.write(aux.to_bytes(1,'little'))
        f = open(self.filename,"a+")
        print(datetime.datetime.now(),file = f)
        print('MessageName:',event.MessageName,file = f)
        print('Key:', event.Key,file = f)
        print('KeyID:', event.KeyID,file = f)
        print('ScanCode:', event.ScanCode,file = f)
        print('---',file = f)
        f.close()


    def shutdown(self):
        self.hm.UnhookKeyboard()
        self.serevent.close() #make sure its closed
        


watcher = Keystroke_Watcher()
pythoncom.PumpMessages()