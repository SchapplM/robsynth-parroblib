% Calculate Gravitation load for parallel robot
% P3RRRRR2G3A0
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% mrSges [4x3]
%   first moment of all robot links (mass times center of mass in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
%
% Output:
% taugX [3x1]
%   forces required to compensate gravitation load
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:13
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RRRRR2G3A0_gravload_para_pf_slag_vp2(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR2G3A0_gravload_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR2G3A0_gravload_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3RRRRR2G3A0_gravload_para_pf_slag_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRRRR2G3A0_gravload_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRRRR2G3A0_gravload_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRRRR2G3A0_gravload_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR2G3A0_gravload_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR2G3A0_gravload_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:12:01
% EndTime: 2020-03-09 21:12:01
% DurationCPUTime: 0.70s
% Computational Cost: add. (609->141), mult. (954->240), div. (72->11), fcn. (612->42), ass. (0->119)
t991 = cos(qJ(3,1));
t1048 = t991 ^ 2;
t988 = cos(qJ(3,2));
t1047 = t988 ^ 2;
t985 = cos(qJ(3,3));
t1046 = t985 ^ 2;
t973 = legFrame(3,2);
t950 = sin(t973);
t953 = cos(t973);
t926 = t953 * g(1) - t950 * g(2);
t972 = mrSges(3,3) - mrSges(2,2);
t949 = g(3) * t972;
t1045 = mrSges(2,1) * t926 + t949;
t974 = legFrame(2,2);
t951 = sin(t974);
t954 = cos(t974);
t927 = t954 * g(1) - t951 * g(2);
t1044 = mrSges(2,1) * t927 + t949;
t975 = legFrame(1,2);
t952 = sin(t975);
t955 = cos(t975);
t928 = t955 * g(1) - t952 * g(2);
t1043 = mrSges(2,1) * t928 + t949;
t996 = mrSges(2,1) * g(3);
t1042 = -t972 * t928 + t996;
t1041 = -t972 * t927 + t996;
t1040 = -t972 * t926 + t996;
t1039 = g(3) * mrSges(1,2);
t987 = cos(qJ(1,3));
t1037 = pkin(1) * t987;
t990 = cos(qJ(1,2));
t1036 = pkin(1) * t990;
t993 = cos(qJ(1,1));
t1035 = pkin(1) * t993;
t1034 = mrSges(3,2) * t926;
t1033 = mrSges(3,2) * t927;
t1032 = mrSges(3,2) * t928;
t977 = sin(qJ(2,3));
t978 = sin(qJ(1,3));
t986 = cos(qJ(2,3));
t923 = t978 * t977 - t987 * t986;
t1031 = t923 * t985;
t980 = sin(qJ(2,2));
t981 = sin(qJ(1,2));
t989 = cos(qJ(2,2));
t924 = t981 * t980 - t990 * t989;
t1030 = t924 * t988;
t983 = sin(qJ(2,1));
t984 = sin(qJ(1,1));
t992 = cos(qJ(2,1));
t925 = t984 * t983 - t993 * t992;
t1029 = t925 * t991;
t976 = sin(qJ(3,3));
t1028 = t950 * t976;
t979 = sin(qJ(3,2));
t1027 = t951 * t979;
t982 = sin(qJ(3,1));
t1026 = t952 * t982;
t1025 = t953 * t976;
t1024 = t954 * t979;
t1023 = t955 * t982;
t957 = 0.1e1 / t977;
t961 = 0.1e1 / t985;
t1022 = t957 * t961;
t958 = 0.1e1 / t980;
t964 = 0.1e1 / t988;
t1021 = t958 * t964;
t959 = 0.1e1 / t983;
t967 = 0.1e1 / t991;
t1020 = t959 * t967;
t1001 = t985 * mrSges(3,1) - t976 * mrSges(3,2);
t969 = qJ(1,3) + qJ(2,3);
t937 = sin(t969);
t940 = cos(t969);
t1019 = t961 * (-(t950 * g(1) + t953 * g(2)) * t1001 + (-t937 * g(3) + t926 * t940) * (mrSges(3,1) * t976 + mrSges(3,2) * t985));
t1000 = t988 * mrSges(3,1) - t979 * mrSges(3,2);
t970 = qJ(1,2) + qJ(2,2);
t938 = sin(t970);
t941 = cos(t970);
t1018 = t964 * (-(t951 * g(1) + t954 * g(2)) * t1000 + (-t938 * g(3) + t927 * t941) * (mrSges(3,1) * t979 + mrSges(3,2) * t988));
t971 = qJ(1,1) + qJ(2,1);
t939 = sin(t971);
t942 = cos(t971);
t999 = t991 * mrSges(3,1) - t982 * mrSges(3,2);
t1017 = t967 * (-(t952 * g(1) + t955 * g(2)) * t999 + (-t939 * g(3) + t928 * t942) * (mrSges(3,1) * t982 + mrSges(3,2) * t991));
t1013 = pkin(2) * t923 * t1046;
t1012 = pkin(2) * t924 * t1047;
t1011 = pkin(2) * t925 * t1048;
t1010 = t986 * t976 * pkin(1);
t1009 = t989 * t979 * pkin(1);
t1008 = t992 * t982 * pkin(1);
t917 = mrSges(3,1) * t926;
t936 = (m(2) + m(3)) * pkin(1) + mrSges(1,1);
t932 = t936 * g(3);
t943 = qJ(3,3) + t969;
t944 = -qJ(3,3) + t969;
t994 = mrSges(3,2) * g(3);
t995 = mrSges(3,1) * g(3);
t908 = (t995 - t1034) * cos(t944) / 0.2e1 + (t917 + t994) * sin(t944) / 0.2e1 + (t995 + t1034) * cos(t943) / 0.2e1 + (t917 - t994) * sin(t943) / 0.2e1 + (t936 * t926 - t1039) * t978 + t1045 * t937 + (mrSges(1,2) * t926 + t932) * t987 + t1040 * t940;
t1007 = t908 * t1022;
t918 = mrSges(3,1) * t927;
t945 = qJ(3,2) + t970;
t946 = -qJ(3,2) + t970;
t909 = (t995 - t1033) * cos(t946) / 0.2e1 + (t918 + t994) * sin(t946) / 0.2e1 + (t995 + t1033) * cos(t945) / 0.2e1 + (t918 - t994) * sin(t945) / 0.2e1 + (t936 * t927 - t1039) * t981 + t1044 * t938 + (mrSges(1,2) * t927 + t932) * t990 + t1041 * t941;
t1006 = t909 * t1021;
t919 = mrSges(3,1) * t928;
t947 = qJ(3,1) + t971;
t948 = -qJ(3,1) + t971;
t910 = (t995 - t1032) * cos(t948) / 0.2e1 + (t919 + t994) * sin(t948) / 0.2e1 + (t995 + t1032) * cos(t947) / 0.2e1 + (t919 - t994) * sin(t947) / 0.2e1 + (t936 * t928 - t1039) * t984 + t1043 * t939 + (mrSges(1,2) * t928 + t932) * t993 + t1042 * t942;
t1005 = t910 * t1020;
t911 = (t1001 * g(3) + t1040) * t940 + t937 * (t1001 * t926 + t1045);
t1004 = t911 * t957 / t1046;
t912 = (t1000 * g(3) + t1041) * t941 + t938 * (t1000 * t927 + t1044);
t1003 = t912 * t958 / t1047;
t913 = (t999 * g(3) + t1042) * t942 + t939 * (t999 * t928 + t1043);
t1002 = t913 * t959 / t1048;
t998 = 0.1e1 / pkin(1);
t997 = 0.1e1 / pkin(2);
t1 = [-g(1) * m(4) + ((-t955 * t1029 + t1026) * t1005 + (-t954 * t1030 + t1027) * t1006 + (-t953 * t1031 + t1028) * t1007) * t998 + (t950 * t1019 + t951 * t1018 + t952 * t1017 + ((t955 * t1011 + (-pkin(2) * t1026 - t955 * t1035) * t991 - t952 * t1008) * t1002 + (t954 * t1012 + (-pkin(2) * t1027 - t954 * t1036) * t988 - t951 * t1009) * t1003 + (t953 * t1013 + (-pkin(2) * t1028 - t953 * t1037) * t985 - t950 * t1010) * t1004) * t998) * t997; -g(2) * m(4) + ((t952 * t1029 + t1023) * t1005 + (t951 * t1030 + t1024) * t1006 + (t950 * t1031 + t1025) * t1007) * t998 + (t953 * t1019 + t954 * t1018 + t955 * t1017 + ((-t952 * t1011 + (-pkin(2) * t1023 + t952 * t1035) * t991 - t955 * t1008) * t1002 + (-t951 * t1012 + (-pkin(2) * t1024 + t951 * t1036) * t988 - t954 * t1009) * t1003 + (-t950 * t1013 + (-pkin(2) * t1025 + t950 * t1037) * t985 - t953 * t1010) * t1004) * t998) * t997; -g(3) * m(4) + (-t937 * t957 * t908 - t938 * t958 * t909 - t939 * t959 * t910 + ((pkin(2) * (t993 * t983 + t984 * t992) * t991 + t984 * pkin(1)) * t913 * t1020 + (pkin(2) * (t990 * t980 + t981 * t989) * t988 + t981 * pkin(1)) * t912 * t1021 + (pkin(2) * (t987 * t977 + t978 * t986) * t985 + t978 * pkin(1)) * t911 * t1022) * t997) * t998;];
taugX  = t1;
