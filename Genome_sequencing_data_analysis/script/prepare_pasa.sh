#!/bin/bash

export LD_LIBRARY_PATH=

sqlDir='./script/mysql_sample'
pasaDir='./script/pasa_sample'

if [ $1 -a -d $1 ];then
    while [ 1 ]
    do
        ((port=$RANDOM%10000+50000))
        if [ ! $( netstat -an|grep LISTEN|grep tcp|awk '{print $4}'|awk -F':' '{if($2)print $2}' | grep $port ) ];then
            break
        fi
    done
    # mysql 
    mkdir -p $1/mysql_bin/tmp
    prodir=`pwd`
    cd $1/mysql_bin
    newDir=`pwd`
    cp -R $sqlDir/bin $sqlDir/data $sqlDir/include $sqlDir/lib $sqlDir/pasa.cnf $newDir
    perl -pi -e 's|<PREFIX>|'$newDir'|' $newDir/pasa.cnf
    perl -pi -e 's|<PORT>|'$port'|' $newDir/pasa.cnf
    cd $newDir/bin && ../bin/mysqld_safe --defaults-file=../pasa.cnf &
    sleep 3
    cd $prodir

    # pasa
    mkdir -p $1/pasa_bin
    cd $1/pasa_bin
    ln -s $pasaDir/* .
    rm pasa_conf;mkdir pasa_conf
    ln -s $pasaDir/pasa_conf/* pasa_conf/
    rm pasa_conf/conf.txt scripts
    cp $pasaDir/pasa_conf/conf.txt pasa_conf/
    cp -R $pasaDir/scripts .
    perl -pi -e 's|MYSQLSERVER=172.17.255.253:3308|MYSQLSERVER=127.0.0.1:'$port'|' pasa_conf/conf.txt
    perl -pi -e 's|USER=pasa|USER=root|' pasa_conf/conf.txt
    perl -pi -e 's|PASSWORD=pasa|PASSWORD=123456|' pasa_conf/conf.txt
fi

