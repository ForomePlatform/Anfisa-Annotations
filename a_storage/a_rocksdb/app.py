#  Copyright (c) 2019. Partners HealthCare and other members of
#  Forome Association
#
#  Developed by Sergey Trifonov based on contributions by Joel Krier,
#  Michael Bouzinier, Shamil Sunyaev and other members of Division of
#  Genetics, Brigham and Women's Hospital
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#        http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
#

import sys, json, logging, signal

from .a_storage import AStorage
from .a_array import AArray

#===============================================
def terminateAll(sig, frame):
    AStorageApp.terminate(sig, frame)

#===============================================
class AStorageApp:
    sConfig = None

    @classmethod
    def setup(cls, config, in_container):
        cls.sConfig = config
        cls.sStorage = AStorage(config)
        cls.sArrays = {name: AArray(cls.sStorage, name, descr)
            for name, descr in config["service"]["arrays"].items()}
        cls.sStorage.activate()

        signal.signal(signal.SIGTERM, terminateAll)
        signal.signal(signal.SIGHUP, terminateAll)
        signal.signal(signal.SIGINT, terminateAll)

    @classmethod
    def terminate(cls, sig, frame):
        cls.sStorage.deactivate()
        logging.info("Application terminated")
        sys.exit(0)

    @classmethod
    def getVersionCode(cls):
        return cls.sVersionCode

    @classmethod
    def hasRunOption(cls, name):
        run_options = cls.sConfig.get("run-options")
        return run_options and name in run_options

    @classmethod
    def getRunModes(cls):
        return cls.sConfig.get("run-modes", [])

    @classmethod
    def getOption(cls, name):
        return cls.sConfig.get(name)

    @classmethod
    def request(cls, serv_h, rq_path, rq_args, rq_descr):
        if rq_path == "/get":
            array_h = cls.sArrays.get(rq_args.get("array"))
            if array_h is None:
                return serv_h.makeResponse("Array not found: "
                    + str(rq_args.get("array")), error = 404)
            report = array_h.request(rq_args)
            return serv_h.makeResponse(mode = "json",
                content = json.dumps(report))
        if rq_path == "/meta":
            return serv_h.makeResponse(mode = "json",
                content = json.dumps(cls.sConfig["service"]["meta"]))
        return serv_h.makeResponse("Page not found",
            error = 404)

    @classmethod
    def checkFilePath(cls, fpath):
        return None
