SDF_PATH='./SDF/IDL'
LARE_PATH='./IDL'
!path = SDF_PATH + PATH_SEP(/SEARCH_PATH) + !path
!path = LARE_PATH + PATH_SEP(/SEARCH_PATH) + !path

@./SDF/IDL/Start.pro
.r getenergy
