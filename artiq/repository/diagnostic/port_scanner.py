#!/usr/bin/env python
import socket
import subprocess
import sys
import argparse

def get_argparser():
    parser = argparse.ArgumentParser(
        description="Scan for open ARTIQ ports.")
    parser.add_argument("ip")
    return parser

artiq_ports = {
    1381: "core device (main)",
    3250: "core device (mon/inj)",
    3248: "InfluxDB bridge",
    3249: "contrller manager",
    3250: "master (notifications)",
    3251: "master (control)",
    3252: "PDQ2",
    3253: "LDA",
    3254: "Novatech 409B",
    3255: "Thorlabs T-Cube",
    3256: "NI PXI6733"
    }

def do_scan(remoteServer):
    remoteServerIP  = socket.gethostbyname(remoteServer)
    try:
        for port in artiq_ports.keys():
            sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            sock.settimeout(0.1)
            result = sock.connect_ex((remoteServerIP, port))
            if result == 0:
                print("({}:{})".format(port, artiq_ports[port]), end="", flush=True)
            else:
                print(".", end="", flush=True)
            sock.close()
        print("")

    except socket.gaierror:
        print('Hostname could not be resolved. Exiting')

    except socket.error:
        print("Couldn't connect to server")

def main():
    args = get_argparser().parse_args()
    do_scan(args.ip)

if __name__ == "__main__":
    main()