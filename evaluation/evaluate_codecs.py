" コーデック評価 "
import glob
import os
import csv
from pathlib import Path
from timeit import timeit
import numpy as np
import scipy.io.wavfile as wv

def _get_wavfile_length_sec(filename):
    """ 音声ファイル長計測 """
    rate, data = wv.read(filename)
    return data.shape[0] / rate

def _measure_execution_time(command):
    """ 実行時間の計測 """
    return timeit(stmt = f'subprocess.run(\'{command}\','\
        'shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)',
        setup = "import subprocess", number = 1)

class Codec:
    """ コーデックインターフェース """
    def __init__(self, compress_option):
        self.compress_option = compress_option

    def get_label(self):
        """ 実験ラベル取得 """
        raise NotImplementedError

    def generate_encode_command(self, in_filename, out_filename):
        """ エンコードコマンド取得 """
        raise NotImplementedError

    def post_encode(self, in_filename, out_filename):
        """ エンコード後の処理（リネームなど） """
        pass

    def generate_decode_command(self, in_filename, out_filename):
        """ デコードコマンド取得 """
        raise NotImplementedError

    def encode(self, in_filename, out_filename):
        """ エンコード """
        escaped_in_filename = f'\"{in_filename}\"'
        command = self.generate_encode_command(escaped_in_filename, out_filename)
        time = _measure_execution_time(command)
        self.post_encode(escaped_in_filename, out_filename)
        return (time, os.path.getsize(out_filename))

    def decode(self, in_filename, out_filename):
        """ デコード """
        command = self.generate_decode_command(in_filename, out_filename)
        return _measure_execution_time(command)

class FLAC(Codec):
    """ FLAC """
    def get_label(self):
        return f"FLAC {self.compress_option}"
    def generate_encode_command(self, in_filename, out_filename):
        return f'flac {self.compress_option} -f -s -o {out_filename} {in_filename}'
    def generate_decode_command(self, in_filename, out_filename):
        return f"flac -d -f -s -o {out_filename} {in_filename}"

class WavPack(Codec):
    """ WavPack """
    def get_label(self):
        return f"WavPack {self.compress_option}"
    def generate_encode_command(self, in_filename, out_filename):
        return f"wavpack {self.compress_option} -q -y {in_filename} -o {out_filename}"
    def generate_decode_command(self, in_filename, out_filename):
        return f"wvunpack {in_filename} -q -y -o {out_filename}"

class TTA(Codec):
    """ TTA """
    def get_label(self):
        return "TTA"
    def generate_encode_command(self, in_filename, out_filename):
        return f"tta -e {in_filename} {out_filename}"
    def generate_decode_command(self, in_filename, out_filename):
        return f"tta -d {in_filename} {out_filename}"

class MonkeysAudio(Codec):
    """ Monkey's Audio """
    def get_label(self):
        return f"Monkey's Audio {self.compress_option}"
    def generate_encode_command(self, in_filename, out_filename):
        return f"mac {in_filename} {out_filename} {self.compress_option}"
    def generate_decode_command(self, in_filename, out_filename):
        # .apeにしないと展開に失敗する...
        in_path = Path(in_filename)
        ape_in_filename = in_path.with_suffix('.ape')
        os.replace(in_filename, ape_in_filename)
        return f"mac {ape_in_filename} {out_filename} -d"

class MPEG4ALS(Codec):
    """ MPEG4-ALS """
    def get_label(self):
        return f"MPEG4-ALS {self.compress_option}"
    def generate_encode_command(self, in_filename, out_filename):
        if os.name == 'nt':
            command = f"mp4alsRM23 {self.compress_option} {in_filename} {out_filename}"
        else:
            command = f"wine64 mp4alsRM23.exe {self.compress_option} {in_filename} {out_filename}"
        return command
    def generate_decode_command(self, in_filename, out_filename):
        if os.name == 'nt':
            command = f"mp4alsRM23 -x {in_filename} {out_filename}"
        else:
            command = f"wine64 mp4alsRM23.exe -x {in_filename} {out_filename}"
        return command

class TAK(Codec):
    """ TAK """
    def get_label(self):
        return f"TAK {self.compress_option}"
    def generate_encode_command(self, in_filename, out_filename):
        # 出力の拡張子を.takにする必要がある
        out_path = Path(out_filename)
        tak_out_filename = out_path.with_suffix('.tak')
        return f'Takc -e -overwrite -silent -tn1 {self.compress_option} {in_filename} {tak_out_filename}'
    def post_encode(self, in_filename, out_filename):
        # .takを目的のファイル名にリネーム
        out_path = Path(out_filename)
        tak_out_filename = out_path.with_suffix('.tak')
        os.replace(tak_out_filename, out_filename)
    def generate_decode_command(self, in_filename, out_filename):
        # .takにしないと展開に失敗する...
        in_path = Path(in_filename)
        tak_in_filename = in_path.with_suffix('.tak')
        os.replace(in_filename, tak_in_filename)
        return f'Takc -d -overwrite -silent -tn1 {tak_in_filename} {out_filename}'

class NARU(Codec):
    """ NARU """
    def get_label(self):
        return f"NARU {self.compress_option}"
    def generate_encode_command(self, in_filename, out_filename):
        return f"naru {self.compress_option} -e {in_filename} {out_filename}"
    def generate_decode_command(self, in_filename, out_filename):
        return f"naru {self.compress_option} -d {in_filename} {out_filename}"

class LINNE(Codec):
    """ LINNE """
    def get_label(self):
        return f"LINNE {self.compress_option}"
    def generate_encode_command(self, in_filename, out_filename):
        return f"linne {self.compress_option} -e {in_filename} {out_filename}"
    def generate_decode_command(self, in_filename, out_filename):
        return f"linne {self.compress_option} -d {in_filename} {out_filename}"

class SRLA(Codec):
    """ SRLA """
    def get_label(self):
        return f"SRLA {self.compress_option}"
    def generate_encode_command(self, in_filename, out_filename):
        return f"srla {self.compress_option} -e {in_filename} {out_filename}"
    def generate_decode_command(self, in_filename, out_filename):
        return f"srla {self.compress_option} -d {in_filename} {out_filename}"

TEST_FILE_DIRECTORY_DICT = {
    "classic": [
        './data/RWC Music Database Sub-Working Group - Rwc-Mdb-C-2001-M*/**/*.wav'
        ],
    "genre": [
        './data/RWC Music Database Sub-Working Group - Rwc-Mdb-G-2001-M*/**/*.wav'
        ],
    "jazz": [
        './data/RWC Music Database Sub-Working Group - Rwc-Mdb-J-2001-M*/**/*.wav'
        ],
    "popular": [
        './data/RWC Music Database Sub-Working Group - Rwc-Mdb-P-2001-M*/**/*.wav'
        ],
    "right": [
        './data/RWC Music Database Sub-Working Group - Rwc-Mdb-R-2001-M*/**/*.wav'
        ],
    }

CODEC_CONFUGURES = [
        # FLAC("-0"),
        # FLAC("-5"),
        # FLAC("-8"),
        # WavPack(""),
        # WavPack("-h"),
        # WavPack("-hh"),
        # WavPack("-x4"),
        # WavPack("-h -x4"),
        # WavPack("-hh -x4"),
        # TTA(""),
        # MonkeysAudio("-c1000"),
        # MonkeysAudio("-c2000"),
        # MonkeysAudio("-c3000"),
        # MonkeysAudio("-c4000"),
        # MPEG4ALS(""),
        # MPEG4ALS("-7"),
        # MPEG4ALS("-b"),
        # LINNE("-m 0"),
        # LINNE("-m 4"),
        # LINNE("-m 7"),
        # NARU("-m 0"),
        # NARU("-m 2"),
        # NARU("-m 4"),
        # TAK("-p0"),
        # TAK("-p0e"),
        # TAK("-p0m"),
        # TAK("-p1"),
        # TAK("-p1e"),
        # TAK("-p1m"),
        # TAK("-p2"),
        # TAK("-p2e"),
        # TAK("-p2m"),
        # TAK("-p3"),
        # TAK("-p3e"),
        # TAK("-p3m"),
        # TAK("-p4"),
        # TAK("-p4e"),
        # TAK("-p4m"),
        SRLA("-m 2 -V 0 -B 4096"),
        SRLA("-m 2 -V 1 -B 4096"),
        SRLA("-m 2 -V 2 -B 4096"),
        SRLA("-m 2 -V 0 -B 8192"),
        SRLA("-m 2 -V 1 -B 8192"),
        SRLA("-m 2 -V 2 -B 8192"),
        SRLA("-m 2 -V 0 -B 16384"),
        SRLA("-m 2 -V 1 -B 16384"),
        SRLA("-m 2 -V 2 -B 16384"),
        SRLA("-m 8 -V 0 -B 4096"),
        SRLA("-m 8 -V 1 -B 4096"),
        SRLA("-m 8 -V 2 -B 4096"),
        SRLA("-m 8 -V 0 -B 8192"),
        SRLA("-m 8 -V 1 -B 8192"),
        SRLA("-m 8 -V 2 -B 8192"),
        SRLA("-m 8 -V 0 -B 16384"),
        SRLA("-m 8 -V 1 -B 16384"),
        SRLA("-m 8 -V 2 -B 16384"),
        SRLA("-m 12 -V 0 -B 4096"),
        SRLA("-m 12 -V 1 -B 4096"),
        SRLA("-m 12 -V 2 -B 4096"),
        SRLA("-m 12 -V 0 -B 8192"),
        SRLA("-m 12 -V 1 -B 8192"),
        SRLA("-m 12 -V 2 -B 8192"),
        SRLA("-m 12 -V 0 -B 16384"),
        SRLA("-m 12 -V 1 -B 16384"),
        SRLA("-m 12 -V 2 -B 16384"),
        SRLA("-m 16 -V 0 -B 4096"),
        SRLA("-m 16 -V 1 -B 4096"),
        SRLA("-m 16 -V 2 -B 4096"),
        SRLA("-m 16 -V 0 -B 8192"),
        SRLA("-m 16 -V 1 -B 8192"),
        SRLA("-m 16 -V 2 -B 8192"),
        SRLA("-m 16 -V 0 -B 16384"),
        SRLA("-m 16 -V 1 -B 16384"),
        SRLA("-m 16 -V 2 -B 16384"),
    ]

if __name__ == "__main__":
    # 一時ファイル名
    COMPRESS_TMP_FILENAME = "compressed.tmp"
    DECOMPRESS_TMP_FILENAME = "decompressd.wav"

    filesdir = {}
    # wavファイルリストを取得
    for c, ds in TEST_FILE_DIRECTORY_DICT.items():
        filesdir[c] = []
        for d in ds:
            filesdir[c] += glob.glob(d, recursive=True)
        # パス文字列を統一的に扱うため'\\'を'/'に置換
        if os.name == 'nt':
            filesdir[c] = [f.replace('\\', '/') for f in filesdir[c]]
        filesdir[c].sort()

    # 計測
    results = {}
    for codec in CODEC_CONFUGURES:
        results[codec.get_label()] = {}
        for categ, fs in filesdir.items():
            results[codec.get_label()][categ]\
                = { 'encode_time': [], 'decode_time': [], 'compress_rate': [] }
            for f in fs:
                print(f'[{codec.get_label()}] {f}')
                # 1ファイル計測
                encode_time, size = codec.encode(f, COMPRESS_TMP_FILENAME)
                decode_time = codec.decode(COMPRESS_TMP_FILENAME, DECOMPRESS_TMP_FILENAME)
                # 結果記録
                original_time = _get_wavfile_length_sec(f)
                results[codec.get_label()][categ]['encode_time']\
                    .append((encode_time * 100) / original_time)
                results[codec.get_label()][categ]['decode_time']\
                    .append((decode_time * 100) / original_time)
                results[codec.get_label()][categ]['compress_rate']\
                    .append((size * 100) / os.path.getsize(f))
                # 中間生成物の削除
                if os.path.exists(COMPRESS_TMP_FILENAME):
                    os.remove(COMPRESS_TMP_FILENAME)
                if os.path.exists(DECOMPRESS_TMP_FILENAME):
                    os.remove(DECOMPRESS_TMP_FILENAME)

    # 全カテゴリ結果を結合した結果の計算
    total_result = {}
    for codec in CODEC_CONFUGURES:
        total_result[codec] = { 'encode_time': [], 'decode_time': [], 'compress_rate': [] }
        for categ in filesdir:
            for e in total_result[codec].keys():
                total_result[codec][e].extend(results[codec.get_label()][categ][e])
        for e in total_result[codec].keys():
            total_result[codec][e] = np.mean(total_result[codec][e])

    # 結果出力
    with open('codec_comparison_summery.csv', 'w', encoding='UTF-8') as f:
        writer = csv.writer(f, lineterminator='\n')
        header = ['']
        for codec in CODEC_CONFUGURES:
            header.append(codec.get_label())
        writer.writerow(header)
        for categ in filesdir:
            row = [f'{categ} mean encode time']
            for codec in CODEC_CONFUGURES:
                row.append(np.mean(results[codec.get_label()][categ]['encode_time']))
            writer.writerow(row)
        row = ['total mean encode time']
        for codec in CODEC_CONFUGURES:
            row.append(total_result[codec]['encode_time'])
        writer.writerow(row)
        for categ in filesdir:
            row = [f'{categ} mean decode time']
            for codec in CODEC_CONFUGURES:
                row.append(np.mean(results[codec.get_label()][categ]['decode_time']))
            writer.writerow(row)
        row = ['total mean decode time']
        for codec in CODEC_CONFUGURES:
            row.append(total_result[codec]['decode_time'])
        writer.writerow(row)
        for categ in filesdir:
            row = [f'{categ} mean compression rate']
            for codec in CODEC_CONFUGURES:
                row.append(np.mean(results[codec.get_label()][categ]['compress_rate']))
            writer.writerow(row)
        row = ['total mean compression rate']
        for codec in CODEC_CONFUGURES:
            row.append(total_result[codec]['compress_rate'])
        writer.writerow(row)
