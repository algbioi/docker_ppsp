#!/usr/bin/env python

"""
    Copyright (C) 2014  Ivan Gregor

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""

import subprocess
import threading
import multiprocessing
import time


class TaskCmd():
    def __init__(self, cmd, cwd='.'):
        """
            Defines one task to be executed on a command line.

            @param cmd: command to be executed on a command line
            @param cwd: current working directory in which the task will be executed
        """
        self.cmd = cmd
        self.cwd = cwd


class _ProcHandlerCmd():
    """
        Stores context of one process.
    """
    def __init__(self, task, timeout=None):
        """
            @type task: TaskCmd
        """
        self.process = subprocess.Popen('exec ' + task.cmd, shell=True, bufsize=-1, cwd=task.cwd)
        self.cmd = task.cmd
        self.runtime = 0.
        self.timeout = timeout

    def incRuntime(self, timeStep):
        self.runtime += timeStep

    def isTimeOut(self):
        if self.timeout is not None and self.runtime > self.timeout:
            return True
        return False

    def getPid(self):
        return self.process.pid


def runCmdParallel(cmdTaskList, maxProc=2, timeout=None, timeStep=0.3):
    """
        Run several command line commands in parallel.

        @param cmdTaskList: list of command line tasks
        @type cmdTaskList: list of TaskCmd
        @param maxProc: maximum number of tasks that will be run in parallel at the same time
        @param timeout: after this number of seconds, the process will be killed, (None if no timeout set)
        @param timeStep: time interval in which processes will be actively checked whether they are running
        @return: list of failed commands, tuple (command, return code, runtime, pid)
    """
    counter = 0
    failList = []
    cmdInGen = True
    cmdGen = iter(cmdTaskList)  # generator of commands
    procArray = {}
    for i in range(maxProc):
        procArray[i] = None

    # loop until all processes finish (or are killed)
    while True:

        # run commands
        if cmdInGen:
            for i in range(maxProc):
                if procArray[i] is None:
                    try:
                        task = cmdGen.next()
                        procArray[i] = _ProcHandlerCmd(task, timeout)  # run process
                        counter += 1
                       # print('Running "%s" cmd: %s' % (counter, task.cmd))
                    except StopIteration:
                        cmdInGen = False  # there are no processes to be run

        # sleep for a while
        time.sleep(timeStep)

        # check for finished processes and processes passed timeout
        for i in range(maxProc):
            ph = procArray[i]
            if ph is not None:  # there is a process in the slot
                if ph.process.poll() is None:
                    # process running
                    ph.incRuntime(timeStep)
                    if ph.isTimeOut():
                        # process over time, kill it!
                        ph.process.kill()
                        print("Process (%s): %s killed! (after %ss)" % (ph.getPid(), ph.cmd, ph.runtime))
                        failList.append((ph.cmd, 9, ph.runtime, ph.getPid()))
                        procArray[i] = None  # free slot
                else:
                    # process finished
                    ph.process.wait()
                    if ph.process.returncode != 0:
                       # print('Process(%s): "%s" ended with return code: "%s' % (ph.getPid(), ph.cmd, ph.process.returncode))
                        failList.append((ph.cmd, ph.process.returncode, ph.runtime, ph.getPid()))
                    procArray[i] = None  # free slot

        # finish if no process is running and there is not process to be run
        if len(set(procArray.values())) == 1 and not cmdInGen:
            break

    return failList


def runCmdSerial(cmdTaskList):
    """
        Run several command line commands one by one

        @param cmdTaskList: list of command line tasks
        @type cmdTaskList: list of TaskCmd
    """
    counter = 0
    for task in cmdTaskList:
        counter += 1
        process = subprocess.Popen('exec ' + task.cmd, shell=True, bufsize=-1, cwd=task.cwd)
        #print('run "%s" cmd: %s' % (counter, task.cmd))
        process.wait()
        #if process.returncode != 0:
           # print('Process: "%s" ended with return code: "%s' % (task.cmd, process.returncode))


class TaskThread():
    def __init__(self, fun, args):
        """
            Defines one function and its arguments to be executed in one thread.

            @param fun: a function to be executed
            @type fun: function
            @param args: arguments of the function
            @type args: tuple
        """
        self.fun = fun
        self.args = args


def runThreadParallel(threadTaskList, maxThreads=2, timeStep=0.05):
    """
        Execute several functions (threads) in parallel.

        @param threadTaskList: list of tasks to be executed
        @type threadTaskList: list of TaskThread
        @param maxThreads: maximum number of tasks that will be run in parallel at the same time
        @param timeStep: time interval in which threads will be actively checked whether they are running
    """
    taskGen = iter(threadTaskList)  # generator of tasks
    taskInGen = True

    threadArray = {}
    for i in range(maxThreads):
        threadArray[i] = None

    # loop until all threads finish
    while True:

        # run commands
        if taskInGen:
            for i in range(maxThreads):
                if threadArray[i] is None:
                    try:
                        task = taskGen.next()
                        t = multiprocessing.Process(target=task.fun, args=task.args)  # init thread
                        t.start()
                        threadArray[i] = t
                    except StopIteration:
                        taskInGen = False  # there are no more threads to be run

        # sleep for a while
        time.sleep(timeStep)

        # check for finished threads
        threadRunning = False
        for i in range(maxThreads):
            t = threadArray[i]
            if t is not None:  # there is a thread in the slot
                if t.is_alive():
                    threadRunning = True
                else:
                    # process finished
                    t.join()
                    threadArray[i] = None  # free slot

        # finish if no thread is running and there is no thread to be run
        if not threadRunning and not taskInGen:
            break


# def f(a, t, n, c):
#     for i in range(n):
#         print(a)
#         cc = 333
#         for j in range(int(10000000*t)):
#             c /= cc
#
#
#
#
# if __name__ == "__main__":
#     f('thread 0', 0.5, 5, 48749394857384234987)
    # runThreadParallel([TaskThread(f, ('thread 1', 0.5, 5, 48749394857384234987)),
    #                     TaskThread(f, ('thread 2', 0.7, 6, 57395769304867332349)),
    #                    TaskThread(f, ('thread 3', 0.8, 7, 87263485768798234987)),
    #                    TaskThread(f, ('thread 4', 0.9, 8, 38573947573957684485)),
    #                    TaskThread(f, ('thread 5', 0.9, 8, 38573947573957684485)),
    #                    TaskThread(f, ('thread 6', 1.0, 8, 38573947573957684485)),
    #                    TaskThread(f, ('thread 7', 1.1, 8, 38573947573957684485))])
