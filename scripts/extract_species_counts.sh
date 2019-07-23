#!/bin/bash

jq '.tree | .. | .data? // empty | .[2]' $1
