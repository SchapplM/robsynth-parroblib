% Calculate Gravitation load for parallel robot
% P3PRRRR8V2G3A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% qJ [3x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
%
% Output:
% taugX [3x1]
%   forces required to compensate gravitation load
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 18:05
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3PRRRR8V2G3A0_gravload_para_pf_slag_vp1(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(8,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V2G3A0_gravload_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V2G3A0_gravload_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P3PRRRR8V2G3A0_gravload_para_pf_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR8V2G3A0_gravload_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRRR8V2G3A0_gravload_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3PRRRR8V2G3A0_gravload_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V2G3A0_gravload_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V2G3A0_gravload_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:04:21
% EndTime: 2020-08-06 18:04:23
% DurationCPUTime: 1.32s
% Computational Cost: add. (927->170), mult. (1986->331), div. (36->4), fcn. (1527->22), ass. (0->124)
t963 = sin(pkin(8));
t1051 = g(3) * t963;
t967 = legFrame(3,2);
t953 = sin(t967);
t956 = cos(t967);
t941 = t956 * g(1) - t953 * g(2);
t965 = cos(pkin(8));
t1002 = -t941 * t965 + t1051;
t971 = sin(qJ(2,3));
t977 = cos(qJ(2,3));
t952 = g(3) * t965;
t1062 = t941 * t963 + t952;
t938 = t953 * g(1) + t956 * g(2);
t964 = sin(pkin(4));
t966 = cos(pkin(4));
t997 = t1062 * t966 - t938 * t964;
t988 = -t1002 * t977 - t997 * t971;
t968 = legFrame(2,2);
t954 = sin(t968);
t957 = cos(t968);
t942 = t957 * g(1) - t954 * g(2);
t1000 = -t942 * t965 + t1051;
t973 = sin(qJ(2,2));
t979 = cos(qJ(2,2));
t1063 = t942 * t963 + t952;
t939 = t954 * g(1) + t957 * g(2);
t996 = t1063 * t966 - t939 * t964;
t987 = -t1000 * t979 - t996 * t973;
t975 = sin(qJ(2,1));
t981 = cos(qJ(2,1));
t969 = legFrame(1,2);
t955 = sin(t969);
t958 = cos(t969);
t943 = t958 * g(1) - t955 * g(2);
t1064 = t943 * t963 + t952;
t940 = t955 * g(1) + t958 * g(2);
t995 = t1064 * t966 - t940 * t964;
t998 = -t943 * t965 + t1051;
t986 = -t995 * t975 - t998 * t981;
t1061 = m(3) / pkin(3);
t970 = sin(qJ(3,3));
t1060 = pkin(2) * t970;
t972 = sin(qJ(3,2));
t1059 = pkin(2) * t972;
t974 = sin(qJ(3,1));
t1058 = pkin(2) * t974;
t976 = cos(qJ(3,3));
t1057 = pkin(3) * t976 ^ 2;
t978 = cos(qJ(3,2));
t1056 = pkin(3) * t978 ^ 2;
t980 = cos(qJ(3,1));
t1055 = pkin(3) * t980 ^ 2;
t1054 = pkin(3) * t976;
t1053 = pkin(3) * t978;
t1052 = pkin(3) * t980;
t1031 = t963 * t964;
t1032 = rSges(3,2) * t964 * t952;
t1040 = t938 * t966;
t1019 = t966 * t970;
t1028 = t964 * t971;
t983 = pkin(7) + pkin(6);
t944 = pkin(2) * t971 - t983 * t977;
t923 = pkin(3) * t1019 + t964 * t944;
t914 = 0.1e1 / (pkin(2) * t1019 + t1028 * t1057 + t923 * t976);
t1050 = (((-t1062 * t964 - t1040) * rSges(3,1) + t988 * rSges(3,2)) * t976 + (t1032 + (t941 * t1031 + t1040) * rSges(3,2) + t988 * rSges(3,1)) * t970) * t914;
t1038 = t939 * t966;
t1017 = t966 * t972;
t1026 = t964 * t973;
t945 = pkin(2) * t973 - t983 * t979;
t924 = pkin(3) * t1017 + t964 * t945;
t915 = 0.1e1 / (pkin(2) * t1017 + t1026 * t1056 + t924 * t978);
t1049 = (((-t1063 * t964 - t1038) * rSges(3,1) + t987 * rSges(3,2)) * t978 + (t1032 + (t942 * t1031 + t1038) * rSges(3,2) + t987 * rSges(3,1)) * t972) * t915;
t1036 = t940 * t966;
t1015 = t966 * t974;
t1024 = t964 * t975;
t946 = pkin(2) * t975 - t983 * t981;
t925 = pkin(3) * t1015 + t964 * t946;
t916 = 0.1e1 / (pkin(2) * t1015 + t1024 * t1055 + t925 * t980);
t1048 = (((-t1064 * t964 - t1036) * rSges(3,1) + t986 * rSges(3,2)) * t980 + (t1032 + (t943 * t1031 + t1036) * rSges(3,2) + t986 * rSges(3,1)) * t974) * t916;
t1010 = m(2) * rSges(2,1) + pkin(2) * m(3);
t950 = (-pkin(6) - rSges(3,3)) * m(3) + m(2) * rSges(2,2);
t1047 = (t988 * t950 + (-t1002 * t971 + t997 * t977) * ((rSges(3,1) * t976 - rSges(3,2) * t970) * m(3) + t1010)) * t914;
t1046 = (t987 * t950 + (-t1000 * t973 + t996 * t979) * ((rSges(3,1) * t978 - rSges(3,2) * t972) * m(3) + t1010)) * t915;
t1045 = (t986 * t950 + (-t998 * t975 + t995 * t981) * ((rSges(3,1) * t980 - rSges(3,2) * t974) * m(3) + t1010)) * t916;
t1044 = t914 * t938;
t1043 = t915 * t939;
t1042 = t916 * t940;
t1030 = t963 * t966;
t1029 = t964 * t970;
t1027 = t964 * t972;
t1025 = t964 * t974;
t1023 = t964 * t976;
t1022 = t964 * t978;
t1021 = t964 * t980;
t1020 = t965 * t966;
t1018 = t966 * t971;
t1016 = t966 * t973;
t1014 = t966 * t975;
t1013 = t966 * t977;
t1012 = t966 * t979;
t1011 = t966 * t981;
t947 = pkin(2) * t977 + t971 * t983;
t1009 = ((-t965 * t1013 + t963 * t971) * t1054 - t947 * t1020 + t944 * t963) * t1050;
t948 = pkin(2) * t979 + t973 * t983;
t1008 = ((-t965 * t1012 + t963 * t973) * t1053 - t948 * t1020 + t945 * t963) * t1049;
t949 = pkin(2) * t981 + t975 * t983;
t1007 = ((-t965 * t1011 + t963 * t975) * t1052 - t949 * t1020 + t946 * t963) * t1048;
t935 = t965 * t1018 + t963 * t977;
t1006 = (t965 * t1023 + t970 * t935) * t1047;
t936 = t965 * t1016 + t963 * t979;
t1005 = (t965 * t1022 + t972 * t936) * t1046;
t937 = t965 * t1014 + t963 * t981;
t1004 = (t965 * t1021 + t974 * t937) * t1045;
t994 = pkin(3) * t1029 - t944 * t966;
t993 = pkin(3) * t1027 - t945 * t966;
t992 = pkin(3) * t1025 - t946 * t966;
t959 = m(1) + m(2) + m(3);
t934 = t963 * t1014 - t965 * t981;
t933 = t963 * t1016 - t965 * t979;
t932 = t963 * t1018 - t965 * t977;
t919 = t965 * t949 + t992 * t963;
t918 = t965 * t948 + t993 * t963;
t917 = t965 * t947 + t994 * t963;
t1 = [-t956 * t1006 - t957 * t1005 - t958 * t1004 - m(4) * g(1) + (-(-(-t955 * t1024 + t934 * t958) * t1055 + (t919 * t958 + t955 * t925) * t980 + (t958 * t1031 + t966 * t955) * t1058) * t1042 - (-(-t954 * t1026 + t933 * t957) * t1056 + (t918 * t957 + t954 * t924) * t978 + (t957 * t1031 + t966 * t954) * t1059) * t1043 - (-(-t953 * t1028 + t932 * t956) * t1057 + (t917 * t956 + t953 * t923) * t976 + (t956 * t1031 + t966 * t953) * t1060) * t1044) * t959 + (t958 * t1007 + t957 * t1008 + t956 * t1009) * t1061; t953 * t1006 + t954 * t1005 + t955 * t1004 - m(4) * g(2) + (-((t958 * t1024 + t934 * t955) * t1055 + (-t919 * t955 + t958 * t925) * t980 + (-t955 * t1031 + t958 * t966) * t1058) * t1042 - ((t957 * t1026 + t933 * t954) * t1056 + (-t918 * t954 + t957 * t924) * t978 + (-t954 * t1031 + t957 * t966) * t1059) * t1043 - ((t956 * t1028 + t932 * t953) * t1057 + (-t917 * t953 + t956 * t923) * t976 + (-t953 * t1031 + t956 * t966) * t1060) * t1044) * t959 + (-t955 * t1007 - t954 * t1008 - t953 * t1009) * t1061; (t963 * t1021 + t974 * t934) * t1045 + (t963 * t1022 + t972 * t933) * t1046 + (t963 * t1023 + t970 * t932) * t1047 - m(4) * g(3) + (-(-t937 * t1055 - t949 * t963 * t980 + (pkin(2) * t1025 + t992 * t980) * t965) * t1042 - (-t936 * t1056 - t948 * t963 * t978 + (pkin(2) * t1027 + t993 * t978) * t965) * t1043 - (-t935 * t1057 - t947 * t963 * t976 + (pkin(2) * t1029 + t994 * t976) * t965) * t1044) * t959 + (((t963 * t1011 + t965 * t975) * t1052 + t949 * t1030 + t946 * t965) * t1048 + ((t963 * t1012 + t965 * t973) * t1053 + t948 * t1030 + t945 * t965) * t1049 + ((t963 * t1013 + t965 * t971) * t1054 + t947 * t1030 + t944 * t965) * t1050) * t1061;];
taugX  = t1;
