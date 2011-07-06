#!/bin/sh
cat patches/*.patch | patch -p0 -b -d wcslib