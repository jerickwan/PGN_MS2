'''
Caching, Counter
env = chem
'''
# %% Globals

import numpy as np
import datetime
import joblib

import os

from collections.abc import Iterable


BOND_CHAR = "."  # used to indicate main chain bonds
BRANCH_CHAR = "--"  # used to indicate side chain bonds /synonyms
shBRANCH_CHAR = "-" # used to indicate side chain bonds /names
BATCH_SIZE = 750  # maximum no. of cpds handled at one time
PRECISION_MASS = 5
PRECISION_INTENSITY = 4
OUTPUT_ADDUCTS = ['[M+H]+', '[M+Na]+', '[M+K]+',
                  '[M+2H]2+', '[M+3H]3+']  # adducts

# %%% Create Folders

for folder in ["img","output","cache"]:
    try:
        cwd = os.getcwd()
        path = cwd+f"/{folder}"
        os.mkdir(path)
    except FileExistsError:
        pass

# %% Developer Options

DEV_MODE = False
if DEV_MODE:
    # enables wakepy, win10toast, telegram
    if __name__ == "__main__":
        print("Wakepy active")
        print("Win10Toast active")
    import wakepy
    from win10toast import ToastNotifier
    TOASTER = ToastNotifier()
    # Telegram
    TELEGRAM_ENABLED = False
    BOT_TOKEN, CHAT_ID, GROUP_ID, TELEBOT = None, None, None, None
    def check_telegram(path="bot.txt"):
        global TELEGRAM_ENABLED, BOT_TOKEN, CHAT_ID, GROUP_ID, TELEBOT
        try:
            import telepot
            BOT_TOKEN, CHAT_ID, GROUP_ID = open(path,"r").read().splitlines()
            TELEBOT = telepot.Bot(BOT_TOKEN)
            TELEGRAM_ENABLED = True
        except Exception:
            TELEGRAM_ENABLED = False
            BOT_TOKEN, CHAT_ID, GROUP_ID, TELEBOT = None, None, None, None
        if __name__ == "__main__" and TELEGRAM_ENABLED:
            print("Telegram active")
        return TELEGRAM_ENABLED
    check_telegram()

else:
    TOASTER = None
    TELEGRAM_ENABLED = False

# %% Common Functions


def make_dir(*lst):
    '''Tries to make a directory for each path in lst.'''
    for dr in lst:
        try:
            os.mkdir(dr)
        except FileExistsError:
            pass
        else:
            print(f"{dr} created.")


def sigmoid(num):
    '''Returns sigmoid of num.'''
    return 1/(1+np.exp(-num))


def relu(num):
    '''Returns ReLu of num'''
    return max(0, num)

def lrelu(num):
    '''Returns Leaky ReLu of num'''
    if num > 0:
        return num
    else:
        return 0.05*num

def flatten(l):
    '''Flattens a multi-layered iterable.'''
    for el in l:
        if isinstance(el, Iterable) and not isinstance(el, (str, bytes)):
            yield from flatten(el)
        else:
            yield el


# %% Caching

CACHE_DIR = "cache"
MEMORY = joblib.Memory(CACHE_DIR, verbose=0,
                       bytes_limit = 5e+9)

# %% Counter

TIME_INIT = datetime.datetime.now()
TIME_STRING = TIME_INIT.strftime("%Y%m%d%H%M")
COUNTER_PRINT_MIN = 250  # min. count to initiate tracking
COUNTER_PRINT_TIME_INTERVAL = 300  # in seconds
COUNTER_PRINT_OBJ_INTERVAL = 500  # in objects


class Counter:

    '''Counts progress of current calculation, outputs progress if necessary.'''

    def __init__(self, min_count=None,
                 time_interval=None, obj_interval=None,
                 broadcast_time=None, chat_id=None):
        """


        Parameters
        ----------
        min_count : integer, optional
            Min. no. of objects to start counter. The default is None.
        time_interval : integer, optional
            Min. no. of seconds elapsed between updates. The default is None.
        obj_interval : integer, optional
            Min. no. of objects counted between updates. The default is None.
        broadcast_time : list, optional
            List of additional events to announce.
            Currently accepts: 'init', 'sleep'
            The default is None.
        chat_id : string, optional
            Chat ID for Telegram messages. If None, defaults to CHAT_ID.

        Returns
        -------
        None.

        """
        self.unit = ""
        self.integer_places = 1
        self.current_time = datetime.datetime.now()
        self.current_count = 0
        self.total_count = 0
        self.initial_time = self.current_time
        # Developer only
        self.wakepy = False # Tracks waking status
        self.toaster = None # Win10Toaster
        self.telegram = None # Telegram
        # Optionals
        if min_count is None:
            self.min_count = COUNTER_PRINT_MIN
        else:
            self.min_count = min_count
        if time_interval is None:
            self.time_interval = COUNTER_PRINT_TIME_INTERVAL
        else:
            self.time_interval = time_interval
        if obj_interval is None or obj_interval == 0:
            self.obj_interval = COUNTER_PRINT_OBJ_INTERVAL
        else:
            self.obj_interval = obj_interval
        if chat_id is None and TELEGRAM_ENABLED:
            self.chat_id = CHAT_ID
        else:
            self.chat_id =  chat_id
        # Broadcast
        if broadcast_time is None or not isinstance(broadcast_time, list):
            self.broadcast_time = []
        else:
            self.broadcast_time = broadcast_time
        if "init" in self.broadcast_time:
            print(f"Init Time : {self.current_time.strftime('%H:%M:%S')}")

    def set_counter(self, total, unit="", current_count=0):
        """
        Starts counting.

        Parameters
        ----------
        total : integer
            Total no. of objects to be counted.
        unit : string, optional
            Name of object. The default is "".
        current_count: integer, optional
            Current count. The default is 0.

        Returns
        -------
        None.

        """
        self.unit = unit
        self.current_count = current_count
        self.current_time = datetime.datetime.now()
        self.initial_time = self.current_time
        self.total_count = total
        self.integer_places = int(np.log10(self.total_count))+1
        self.print_current_count()
        self.disable_sleep()

    def update_boolean(self):
        """
        Checks if an update needs to be outputted.

        Returns
        -------
        boolean
            If True, update to be outputted.

        """
        if (self.total_count > self.min_count or
                self.total_count == 0):
            self.current_time = datetime.datetime.now()
            time_interval = self.current_time - self.initial_time
            return time_interval.seconds > self.time_interval\
                or self.current_count % self.obj_interval == 0
        else:
            return False

    def update_counter(self):
        """
        Update internal count and prints updates intermittently.

        Returns
        -------
        None.

        """
        self.current_count += 1
        if self.update_boolean():
            self.initial_time = self.current_time
            self.print_current_count()
        # Check completion
        if self.current_count == self.total_count:
            self.enable_sleep()

    def print_current_count(self):
        """
        Prints update.

        Returns
        -------
        None.

        """
        time = self.current_time.strftime("%H:%M:%S")
        current_count = f"{self.current_count:0{self.integer_places}d}"
        print(
            f"\t\t{current_count}/{self.total_count} {self.unit}\t...\t{time}")

    def disable_sleep(self):
        """
        Disables sleep.

        Returns
        -------
        None.

        """
        if not DEV_MODE:
            return
        elif not self.wakepy:
            wakepy.set_keepawake(keep_screen_awake=False)
            self.wakepy = True
            if "sleep" in self.broadcast_time:
                print("\tSleeping disabled.")

    def enable_sleep(self):
        """
        Enables sleep.

        Returns
        -------
        None.

        """
        if not DEV_MODE:
            return
        if self.wakepy:
            wakepy.unset_keepawake()
            self.wakepy = False
            if "sleep" in self.broadcast_time:
                print("\tSleeping enabled.")

    def show_toast(self, title, msg, duration, **kwargs):
        """
        Displays a Windows Toast notification.

        Parameters
        ----------
        title : string
        msg : string
        duration : float
            Duration in seconds.
        **kwargs :
            Passed to ToastNotifier.show_toast()

        Returns
        -------
        None.

        """
        if not DEV_MODE:
            return
        if self.toaster is None:
            self.toaster = TOASTER
        try:
            self.toaster.show_toast(title, msg, duration=duration, threaded=True,
                                    **kwargs)
        except Exception as exc:
            print(exc)

    def send_telegram_msg(self, msg, chat_id=None):
        """
        Sends a Telegram message.

        Parameters
        ----------
        msg : string
        chat_id : string, optional
            ID of user. If None, the default is self.chat_id.

        Returns
        -------
        None.

        """
        if not DEV_MODE:
            return
        if TELEGRAM_ENABLED:
            if self.telegram is None:
                self.telegram = TELEBOT
            if chat_id is None:
                chat_id = self.chat_id
            try:
                self.telegram.sendMessage(chat_id,msg)
            except Exception as exc:
                print(exc)
