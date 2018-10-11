// arif.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <stdio.h>
#include <stdlib.h>


struct table_
{
    unsigned char symbol;
    int number;
};

class bit_file
{
    public:
    FILE *_file;
    int current_bit;
    unsigned char _temp;
    unsigned char _mask;
    bool _EOF;
    bit_file (const char* file_name, const char* mode);
    ~bit_file ();
    int ReadBit ();
    bool EoF ();
    void WriteBit (int bit);
    char ReadSymbol ();
    int  ReadByte ();
    void Writesym (char c);
    void WriteSize (int size);
    int  ReadSize ();
    void CompleteByte ();
    int  ReadTable (table_*& table);
    void WriteTable (table_* table, unsigned char size);
    int  CreateTable (table_*& table, int& table_size, int* IndexTable);
    void rewind_ ();
};

void Compress (char* comp_name, char* data_name);
void Decompress (char* comp_name, char* data_name);

int main(int argc, char* argv[])
{
	if (argv[1][0] == 'c')
	{
		Compress (argv [3], argv[2]);
	}
	if (argv[1][0] == 'd')
	{
		Decompress (argv [2], argv [3]);
	}
    return 0;
}

void Compress (char* comp_name, char* data_name)
{
    bit_file DataFile (data_name, "rb");
    bit_file CompressedFile (comp_name, "wb");
    unsigned char c;
    int j;
    unsigned long long int l;
    unsigned long long int li = 0;
    unsigned long long int hi = 0xFFFFFF;
    int IndexTable [256] = {0};
    int table_size;
    table_* table;
    int file_size = DataFile.CreateTable (table, table_size, IndexTable);
    int delitel   = table[table_size - 1].number;
    int First_qtr = (hi +1)/4;
    int Half      = First_qtr*2;
    int Third_qtr = First_qtr*3;
    int bits_to_follow = 0;
    CompressedFile.WriteSize (file_size);
    CompressedFile.WriteTable (table, table_size);
    printf ("%d file_size\n", file_size);
    printf ("%d table_size\n", table_size);
    DataFile.rewind_();
    for (int i = 0; i < file_size; i++)
    {
        c = DataFile.ReadSymbol();
        j = IndexTable[c];
        l = li;
        li = li + (table[j-1].number)*(hi - li + 1)/delitel;
        hi = l + (table[j].number)*(hi - l + 1)/delitel - 1;
        for(;;)
        {
            if(hi < Half)
            {
                CompressedFile.WriteBit(0);
                for(; bits_to_follow > 0; bits_to_follow--)
                    CompressedFile.WriteBit(1);
            }
            else if(li >= Half)
            {
                CompressedFile.WriteBit(1);
                for(; bits_to_follow > 0; bits_to_follow--)
                    CompressedFile.WriteBit(0);
                li -= Half;
                hi -= Half;
            }
            else if((li > First_qtr)&&(hi < Third_qtr))
            {
                bits_to_follow++;
                li -= First_qtr;
                hi -= First_qtr;
            }
            else break;
            li += li;
            hi += hi + 1;
        }
    }
    CompressedFile.WriteBit(1);
    CompressedFile.CompleteByte();
}

void Decompress(char* comp_name, char* data_name)
{
    int j;
    unsigned long long int l;
    unsigned long long int li = 0;
    unsigned long long int hi = 0xFFFFFF;
    char c;
    bit_file CompressedFile (comp_name, "rb");
    bit_file DataFile (data_name, "wb");
    table_* table;
    int file_size  = CompressedFile.ReadSize ();
    int table_size = CompressedFile.ReadTable (table);
    int delitel    = table[table_size - 1].number;
    int First_qtr  = (hi +1)/4;
    int Half       = First_qtr*2;
    int Third_qtr  = First_qtr*3;
    int value      = (CompressedFile.ReadByte())*256;
    value += CompressedFile.ReadByte();
    value <<= 8;
    value += CompressedFile.ReadByte();
    for (int i = 0; i < file_size; i++)
    {
		j = 1;
        while (value > (li + (table[j].number)*(hi - li + 1)/delitel - 1))
        {
            j++;
        }
        c = table[j].symbol;
        l = li;
        li = li + table[j - 1].number*(hi - li + 1)/delitel;
        hi = l + table[j].number*(hi - l + 1)/delitel - 1;
        for(;;)
        {
            if(hi < Half)
					;
            else if(li >= Half)
            {
                value -= Half;
                li -= Half;
                hi -= Half;
            }
            else if((li > First_qtr)&&(hi < Third_qtr))
            {
                value -= First_qtr;
                li -= First_qtr;
                hi -= First_qtr;
            }
            else break;
            li += li;
            hi += hi +1;
            value += value + CompressedFile.ReadBit();
        }
        DataFile.Writesym(c); 
    }
}

int bit_file::CreateTable (table_*& table, int& table_size, int* IndexTable)
{
    int b[256] = {0};
    unsigned char temp;
    int size = 0;
    table_size = 256;
    while (true)
    {
        if (fread (&temp, 1, 1, _file) != 1) break;
        b[temp]++;
        size++;
    }
    for (int i = 0; i < 256; i++)
        if (b[i] == 0) table_size --;
    table = new table_[++table_size];
    int i = 1;
    table[0].number = 0;
    table[0].symbol = 0;
    for (int j = 0; j < 256; j++)
    {
        if (b[j] != 0)
        {
            table[i].number = b[j];
            table[i].symbol = j;
            i++;
        }
    }
    table_ temp_t;
    for (int i = 1; i < table_size - 1; i++)
        for (int j = i + 1; j < table_size; j++)
        {
            if (table[j].number > table[i].number)
            {
                temp_t = table[j];
                table[j] = table[i];
                table[i] = temp_t;
            }
        };

    for (int i = 1; i < table_size; i++)
    {
        IndexTable[table[i].symbol] = i;
    }
    for (int i = 1; i < table_size; i++)
        table[i].number += table[i - 1].number;
    return size;
}

void bit_file::WriteTable (table_* table, unsigned char size)
{
    fwrite (&size, 1, 1, _file);
    for (int i = 1; i < size; i++)
    {
        fwrite (&(table[i].symbol), 1, 1, _file);
        fwrite (&(table[i].number), sizeof(int), 1, _file);
    }
}

int bit_file::ReadTable (table_*& table)
{
   unsigned char size;
    fread (&size, 1, 1, _file);
    table = new table_[size];
    table[0].number = 0;
    table[0].symbol = 0;
    for (int i = 1; i < size; i++)
    {
        fread (&(table[i].symbol), 1, 1, _file);
        fread (&(table[i].number), sizeof(int), 1, _file);
    }
    return size;
}

bit_file::~bit_file ()
{
    fclose (_file);
}

bit_file::bit_file (const char* file_name, const char* mode)
{
    _file = fopen (file_name, mode);
    current_bit = 0;
    _temp = 0;
    _mask = 0x80;
    _EOF = false;
}

int bit_file::ReadBit ()
{
    if (current_bit == 0)
    {
        _temp = this->ReadByte();
        _mask = 0x80;
    }
    current_bit = (current_bit + 1) % 8;
    if ((_temp & _mask) == 0)
    {
        _mask >>= 1;
        return 0;
    }
    else
    {
        _mask >>= 1;
        return 1;
    }
}

bool bit_file::EoF ()
{
    return _EOF;
}

char bit_file::ReadSymbol ()
{
    char temp;
    fread(&temp, 1, 1, _file);
    return temp;
}

int bit_file::ReadByte ()
{
    unsigned char temp;
    if (_EOF) return 0;
    _EOF = fread(&temp, 1, 1, _file) != 1;
    if (_EOF) return 0;
    return temp;
}

void bit_file::WriteBit (int bit)
{
    _temp += _mask*bit;
    current_bit = (current_bit + 1) % 8;
    if (current_bit == 0)
    {
        fwrite (&_temp, 1, 1, _file);
        _temp = 0;
        _mask = 0x80;
    }
    else _mask >>= 1;
}

void bit_file::Writesym (char c)
{
    fwrite (&c, 1, 1, _file);
}

void bit_file::CompleteByte ()
{
    if (current_bit != 0)
    {
        fwrite (&_temp, 1, 1, _file);
    }
}

void bit_file::WriteSize (int size)
{
    fwrite (&size, sizeof (int), 1, _file);
}

int bit_file::ReadSize ()
{
    int size;
    fread (&size, sizeof (int), 1, _file);
    return size;
}

void bit_file::rewind_ ()
{
    rewind (_file);
    _EOF = false;
}