/*
 * Written by Bastien Chevreux (BaCh)
 *
 * Copyright (C) 2000 and later by Bastien Chevreux
 *
 * All rights reserved.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the
 * Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA
 *
 */

#ifndef _timerestrict_h_
#define _timerestrict_h_

#ifdef TIMERESTRICTED
#define TR_OUT_MAXYEAR 2008
#define TR_OUT_MAXMONTH 12
#define TR_OUT_MAXDAY 13

#define TR_MAXYEAR (TR_OUT_MAXYEAR-1900)
#define TR_MAXMONTH (TR_OUT_MAXMONTH-1)
#define TR_MAXDAY 13
#endif

#endif
