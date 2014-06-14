// Tested with Arduino 1.05
// polling RFM12B to decode FSK iT+ with a Jeenode from Jeelabs.
// Arduino UNO + RFM12B Module
// device  
//         : IT+ TX22U-IT          LA Crosse Technology (c) USA
// info    : http://www.g-romahn.de/ws1600/Datepakete_raw.txt
//           http://forum.jeelabs.net/node/110
//           http://fredboboss.free.fr/tx29/tx29_sw.php
//           http://www.f6fbb.org/domo/sensors/
//           http://www.mikrocontroller.net/topic/67273 benedikt.k org
//           rinie,marf,joop,Mohammad 23 Aug 2013

#include <avr/io.h>
#include <avr/interrupt.h>
#include <util/delay.h>

#define clrb(sfr, bit)     (_SFR_BYTE(sfr) &= ~_BV(bit)) 
#define setb(sfr, bit)     (_SFR_BYTE(sfr) |= _BV(bit)) 

#define RF_PORT	PORTB
#define RF_DDR	DDRB
#define RF_PIN	PINB

#define SDI   3              // polling  4 Pins ,3 SPI + Chipselect
#define SCK1  5              // differs with Jee-node lib !
#define CS    2  
#define SDO   4 

// *******************************************************************************
static uint8_t quartets[26];
void setup () 
    {
    Serial.begin(9600);

    rf12_la_init();
    }

void loop () 
          {
            rf12_xfer(0xA534);
            receive();
            rf12_xfer(0xA7d0); 
            receive();
            rf12_xfer(0xAA6b);
            receive();
          }

void receive(void)                         
{	unsigned char test[13];	  
        uint8_t i,temp;
   
	rf12_rxdata(test,13);
             
        for(i = 0; i<13; i++){                                
            temp = test[i];
            quartets[2*i+1]= temp & 0b00001111;
            quartets[2*i]= (temp>>4) & 0b00001111;
           }
  	
           display_values();
           Serial.println(); 
	
}

void rf12_rxdata(unsigned char *data, unsigned char number)
{	uint8_t  i;
	rf12_xfer(0x82C8);			// receiver on
	rf12_xfer(0xCA81);			// set FIFO mode
	rf12_xfer(0xCA83);			// enable FIFO
	for (i=0; i<number; i++)
       
	{	rf12_ready();
      
      
		*data++=rf12_xfer(0xB000);
	}
	rf12_xfer(0x8208);			// Receiver off 
}

// ******** SPI + RFM 12B functions   ******************************************

unsigned short rf12_xfer(unsigned short value)
{	uint8_t i;

	clrb(RF_PORT, CS);
	for (i=0; i<16; i++)
	{	if (value&32768)
			setb(RF_PORT, SDI);
		else
			clrb(RF_PORT, SDI);
		value<<=1;
		if (RF_PIN&(1<<SDO))
			value|=1;
		setb(RF_PORT, SCK1);
		asm("nop");
		asm("nop");
		clrb(RF_PORT, SCK1);
	}
	setb(RF_PORT, CS);
	return value;
}


void rf12_ready(void)
{
clrb(RF_PORT, CS);
asm("nop");
asm("nop");
while (!(RF_PIN&(1<<SDO))); // wait until FIFO ready
setb(RF_PORT, CS);
}


static void rf12_la_init() 
{
RF_DDR=(1<<SDI)|(1<<SCK1)|(1<<CS);
RF_PORT=(1<<CS);
for (uint8_t  i=0; i<10; i++) _delay_ms(10); // wait until POR done
rf12_xfer(0x80f7); // 80e8 CONFIGURATION EL,EF,868 band,12.0pF      // iT+ 915  80f7 
rf12_xfer(0xA534); // a67c FREQUENCY SETTING 909.99//                 
rf12_xfer(0xC627); // c627 DATA RATE  8.621 kbps
rf12_xfer(0xC26a); // c26a DATA FILTER COMMAND 
rf12_xfer(0xCA12); // ca12 FIFO AND RESET  8,SYNC,!ff,DR 
rf12_xfer(0xCEd4); // ced4 SYNCHRON PATTERN  0x2dd4 
rf12_xfer(0xC49f); // c49f AFC during VDI HIGH +15 -15 AFC_control_commAND
rf12_xfer(0x95c0); // 95c0 RECEIVER CONTROL VDI Medium 67khz LNA max DRRSI 103 dbm  
rf12_xfer(0xCC77); // cc77 not in RFM01 
rf12_xfer(0x9872); // 9872 transmitter not checked
rf12_xfer(0xE000); // e000 NOT USE 
rf12_xfer(0xC800); // c800 NOT USE 
rf12_xfer(0xC040); // c040 1.66MHz,2.2V 
}

void display_values(void)
{

int temp;
uint8_t i;

if ( (quartets[0] != 0xa) ) return;	// if bad data no display

for ( i = 1; i < quartets[3]+1; i++) {
	
	switch( quartets[i*4] ) {
	case 0:
		//temperature
		temp = (quartets[i*4+1]*100 + quartets[i*4+2]*10 + quartets[i*4+3])-400;
		Serial.print(temp/10);
                Serial.print(".");
                Serial.print(quartets[i*4+3],HEX);
                Serial.print("\t");
	break;
	case 1:
		 //humidity
                 Serial.print("\t");
		 Serial.print(quartets[i*4+2],HEX);
		 Serial.print(quartets[i*4+3],HEX);
	         Serial.print(" ");
                 Serial.print("\t");
                 Serial.print("\t");
	break;
	case 2:
		 //Rain
		 temp = ( quartets[i*4+3] + (quartets[i*4+2]<<4) + ((int)quartets[i*4+1]<<8) )*0.5;
		 Serial.print(temp);
                 Serial.print("\t");
                 Serial.print("\t");
	break;
	case 3:
		 //Wind
		temp = ( quartets[i*4+3] + (quartets[i*4+2]<<4) );
		 Serial.print(temp);
                 Serial.print("\t");
                 Serial.print("\t");
                 
		 //Wind Direction
		temp = (quartets[i*4+1]*225/10);
		 Serial.print(temp);
                 Serial.print("\t");
                 Serial.print("\t");
	break; 
	case 4:
		 //Gust
		temp = ( quartets[i*4+3] + (quartets[i*4+2]<<4) + ((int)quartets[i*4+1]<<8) )*36/10;
		Serial.print(temp);
                
	break; 
	}    
}   
}


