;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; read_ncdf.pro - This file contains IDL functions to read netCDF data files
;               into IDL variables.
;
;  This file contains the following functions and procedures:
;     functions:
;        getFile - strips the directory and optionally) any suffixes from
;                  the path+file
;        getDir - returns the directory from the full path+file
;        validateName - validates the name that will be used as a netCDF
;                       variable
;     procedures:
;        read_ncdf - reads netCDF variables into IDL
;
;  History:
;  Date       Name          Action
;  ---------  ------------  ----------------------------------------------
;  06 Jun 97  S. Rupert     Created.
;  09 Jun 97  S. Rupert     Fully tested.
;  10 Jun 97  S. Rupert     Modified keyword usage.
;  03 Feb 98  S. Rupert     Added additional error checking, and warning to
;                           output script.
;  17 Feb 98  S. Rupert     Corrected validation routine to handle instance
;                           of name strating with a number and containing a
;                           dash.
;  17 Mar 02  D. Morrison   Changed routine to return data into a structure
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


function getFile, fullpath, suffix=suffix
; func_description
; This function returns the filename name from the full path.
; Inputs:    fullpath - full directory+file path
; Keyword:   /suffix:   include inptu suffix in output file name
; Outputs:   file - filename
; Example Call: file = getFile(fullpath)

; Retrieve the postion at which the first '/' character occurs from
; the end of the string.
dirlen = strpos(fullpath, '/')

; Retrieve the full length of the original string.
len = strlen(fullpath)

; Retrieve the filename.
fullfile = strmid(fullpath, dirlen+1, len)

; Retrieve the position at which the first '.' character occurs from
; the end of the string.
len = -1
if not(keyword_set(suffix)) then len = strpos(fullfile, '.')
if (len EQ -1) then len = strlen(fullfile)

; Retrieve the file.
file = strmid(fullfile, 0, len)

; Return the file name.
return, file

; End function.
end


function getDir, fullpath
; func_description
; This function returns the directory name from the full path.
; Inputs:   fullpath - full directory+file path
; Outputs:  dir - directory path
; Example Call: dir = getDir(fullpath)

; Retrieve the postion at which the first '/' character occurs from
; the end of the string.
len = strpos(fullpath, '/')

; Retrieve the filename.
if (len EQ -1) then dir = "./" $
else dir = strmid(fullpath, 0, len+1)

; Return the file name.
return, dir

; End function.
end


function validateName, varname
; func_description
; This routine ensures that the given name does not start with a number,
; nor contain a dash.  IDL cannot accept a variable starting with a
; number or containing a dash.  If the name starts with a number, an
; underscore is prepended to the name, and if it contains a dash, the
; dash is replaced with an underscore.

; Initialize the name.
name = varname

; If the name starts with a number, prepend it with an underscore.
if (strpos(varname, '0') EQ 0) then name = strcompress("_"+varname)
if (strpos(varname, '1') EQ 0) then name = strcompress("_"+varname)
if (strpos(varname, '2') EQ 0) then name = strcompress("_"+varname)
if (strpos(varname, '3') EQ 0) then name = strcompress("_"+varname)
if (strpos(varname, '4') EQ 0) then name = strcompress("_"+varname)
if (strpos(varname, '5') EQ 0) then name = strcompress("_"+varname)
if (strpos(varname, '6') EQ 0) then name = strcompress("_"+varname)
if (strpos(varname, '7') EQ 0) then name = strcompress("_"+varname)
if (strpos(varname, '8') EQ 0) then name = strcompress("_"+varname)
if (strpos(varname, '9') EQ 0) then name = strcompress("_"+varname)

; If the name contains a dash replace it with an underscore.
if (strpos(name, '-') NE -1) then begin
   pieces = str_sep(name, '-')
   n_pieces = n_elements(pieces)
   name = pieces(0)
   for i=1,n_pieces-1 do begin
      name = strcompress(name+"_"+pieces(i))
   endfor
endif

; Return the file name.
return, name

; End function.
end


pro read_ncdf, infile, data, attrib=attrib
; pro_description, infile
; This procedure reads the data in a given netCDFfile into an IDL data structure
; Inputs:           infile  - full path to netCDF file of interest
;                   attrib - includes extractions of all input file
;                             attributes in idl script
; Outputs:          data - structure into which data is returned
; Example:  read_ncdf, '/vol1/cids/netCDF.file', data, /attrib

; Establish that the correct data have been input.
if n_params() LT 1 then begin
   print, "Usage:  read_ncdf, 'inputfile', data, "+$
          "/attrib"+string(10B)+$
          string(9B)+"inputfile  - full path to netCDF file of interest"+string(10B)+$
          string(9B)+"data - desired name of resultant idl data structure"+string(10B)+$
          string(9B)+"/attrib:    keyword indicating inclusion of netCDF attributes"
   return
endif

; Ensure that the netCDF format is supported on the current platform.
if not(ncdf_exists()) then begin
   print, "The Network Common Data Format is not supported on this platform."
   return
endif

; Open the netcdf file for reading.
ncid = NCDF_OPEN(infile)
if (ncid EQ -1) then begin
   print, "The file "+infile+" could not be opened, please check the path."
   return
endif

; Retrieve general information about this netCDF file.
ncidinfo = NCDF_INQUIRE(ncid)

; Write the file header.
dir = getDir(infile)
if (dir EQ './') then dir = 'the current directory'
file = getFile(infile,/suffix)

; Place the desired variables in local arrays.
for i=0, ncidinfo.Nvars-1 do begin
   vardata = NCDF_VARINQ(ncid, i)
   valid_varname = validateName(vardata.Name)
   varname=valid_varname
   NCDF_VARGET, ncid, i, varname      ; Read in variable vardata.Name
   if i eq 0 then begin
 		data = create_struct(valid_varname,varname)
   		endif else begin
   		  data = create_struct(data, valid_varname,varname)
   		endelse
   if (keyword_set(attrib)) then begin
      for j=0, vardata.Natts-1 do begin
         att = NCDF_ATTNAME(ncid, i, j)
         attname = strcompress(valid_varname+"_"+strcompress(att,/REMOVE_ALL))
	print,attname
         valid_attname = attname
         NCDF_ATTGET, ncid, i,att,attname
         attname = STRING(attname)
   		 data = create_struct(data, valid_attname,attname)
      endfor
   endif
endfor

if (keyword_set(attrib)) then begin
  ; Read in the Global Attributes."
  for i=0, ncidinfo.Ngatts-1 do begin
     name = NCDF_ATTNAME(ncid, /GLOBAL, i)
     attname = validateName(name)
	 valid_attname = attname
     NCDF_ATTGET, ncid, /GLOBAL, name, attname
     attname = STRING(attname)
     data = create_struct(data, valid_attname, attname)
  endfor
 endif

NCDF_CLOSE, ncid      ; Close the NetCDF file"

; Return to the caller.
return

; End procedure.
end
