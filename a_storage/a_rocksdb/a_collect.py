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
import json
from pyliftover import LiftOver
#===============================================
class ACollect:
    def __init__(self, arrays):
        self.mSchemaSeq_hg19 = []
        self.mSchemaSeq_hg38 = []
        self.mFastaSchema = None
        self.mFastaTypes = None
        for array_h in arrays.values():
            if array_h.getDBKeyType() == "hg19":
                for schema_h in array_h.iterSchemaSeq():
                    self.mSchemaSeq_hg19.append(schema_h)
            elif array_h.getDBKeyType() == "hg38":
                for schema_h in array_h.iterSchemaSeq():
                    self.mSchemaSeq_hg38.append(schema_h)
            else:
                assert array_h.getDBKeyType() == "fasta"
                for schema_h in array_h.iterSchemaSeq():
                    assert self.mFastaSchema is None
                    self.mFastaSchema = schema_h
        self.mSchemaNames_hg19 = {schema_h.getName()
            for schema_h in self.mSchemaSeq_hg19}
        self.mSchemaNames_hg38 = {schema_h.getName()
            for schema_h in self.mSchemaSeq_hg38}
        self.mAllSchemaNames = self.mSchemaNames_hg19 | self.mSchemaNames_hg38
        self.mLO_19_38 = LiftOver('hg19', 'hg38')
        self.mLO_38_19 = LiftOver('hg38', 'hg19')
        if self.mFastaSchema is not None:
            self.mFastaTypes = {tp_name: "fasta/" + tp_name
                for tp_name in ("hg19", "hg38")}
            assert (sorted(self.mFastaSchema.iterTypeNames())
                == sorted(self.mFastaTypes.keys()))
            self.mAllSchemaNames |= set(self.mFastaTypes.values())

    def _conv19_38(self, chrom, pos):
        coord = self.mLO_19_38.convert_coordinate(chrom, pos - 1)
        if not coord:
            return None
        return coord[0][1] + 1

    def _conv38_19(self, chrom, pos):
        coord = self.mLO_38_19.convert_coordinate(chrom, pos - 1)
        if not coord:
            return None
        return coord[0][1] + 1

    def request(self, rq_args, rq_descr):
        if "@request" in rq_args:
            rq_info = rq_args["@request"]
            variants = rq_info["variants"]
            fasta = rq_info.get("fasta", "hg38")
            arrays = rq_info.get("arrays", self.mAllSchemaNames)
        else:
            variants = json.loads(rq_args["variants"])
            fasta = rq_args.get("fasta", "hg38")
            if "arrays" in rq_args:
                arrays = json.loads(rq_args["arrays"])
            else:
                arrays = self.mAllSchemaNames
        if not isinstance(arrays, set):
            arrays = set(arrays)
        result = []
        for var_info in variants:
            chrom, pos = var_info["chrom"], var_info["pos"]
            ret = {"chrom": chrom, "pos": pos}
            last_pos = None
            if self.mFastaSchema is not None:
                last_pos = var_info.get("last")
            if last_pos is not None:
                ret["last"] = last_pos
            result.append(ret)
            filtering = dict()
            for key in ("alt", "ref"):
                val = var_info.get(key)
                if val:
                    ret[key] = val
                    filtering[key] = val
            if len(filtering) == 0:
                filtering = None
            if len(arrays & self.mSchemaNames_hg19) > 0:
                if fasta == "hg38":
                    _pos = self._conv38_19(chrom, pos)
                    ret["hg19"] = _pos
                    if last_pos is not None:
                        _last_pos = self._conv38_19(chrom, last_pos)
                        ret["hg19-last"] = _last_pos
                else:
                    _pos, _last_pos = pos, last_pos

                if _pos is not None:
                    key = (chrom, _pos)
                    for schema_h in self.mSchemaSeq_hg19:
                        if schema_h.getName() not in arrays:
                            continue
                        rec_data = schema_h.getRecord(key, filtering)
                        if rec_data is not None:
                            ret[schema_h.getName()] = rec_data
                    if (last_pos is not None
                            and self.mFastaTypes["hg19"] in arrays):
                        ret[self.mFastaTypes["hg19"]] = (
                            self.mFastaSchema.getRecord(
                                key, {"type": "hg19"}, _last_pos))

            if len(arrays & self.mSchemaNames_hg38) > 0:
                if fasta == "hg19":
                    _pos = self._conv19_38(chrom, pos)
                    ret["hg38"] = _pos
                    if last_pos is not None:
                        _last_pos = self._conv19_38(chrom, last_pos)
                        ret["hg38-last"] = _last_pos
                else:
                    _pos, _last_pos = pos, last_pos
                if _pos is not None:
                    key = (chrom, _pos)
                    for schema_h in self.mSchemaSeq_hg38:
                        if schema_h.getName() not in arrays:
                            continue
                        rec_data = schema_h.getRecord(key, filtering)
                        if rec_data is not None:
                            ret[schema_h.getName()] = rec_data
                        if (last_pos is not None
                                and self.mFastaTypes["hg38"] in arrays):
                            ret[self.mFastaTypes["hg38"]] = (
                                self.mFastaSchema.getRecord(
                                    key, {"type": "hg38"}, _last_pos))
        return result
