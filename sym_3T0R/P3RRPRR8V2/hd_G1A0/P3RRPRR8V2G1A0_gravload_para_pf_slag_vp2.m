% Calculate Gravitation load for parallel robot
% P3RRPRR8V2G1A0
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
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
% Datum: 2022-11-07 13:12
% Revision: e482436b586c4f286726c907c195760c5ac72455 (2022-11-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RRPRR8V2G1A0_gravload_para_pf_slag_vp2(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR8V2G1A0_gravload_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR8V2G1A0_gravload_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RRPRR8V2G1A0_gravload_para_pf_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR8V2G1A0_gravload_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR8V2G1A0_gravload_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRPRR8V2G1A0_gravload_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR8V2G1A0_gravload_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR8V2G1A0_gravload_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-11-07 13:11:21
% EndTime: 2022-11-07 13:11:21
% DurationCPUTime: 0.42s
% Computational Cost: add. (630->124), mult. (783->191), div. (27->6), fcn. (543->35), ass. (0->106)
t1034 = qJ(3,1) + pkin(5);
t1033 = qJ(3,2) + pkin(5);
t1032 = qJ(3,3) + pkin(5);
t1025 = qJ(2,3) + pkin(7);
t1081 = pkin(3) * cos(t1025);
t1026 = qJ(2,2) + pkin(7);
t1080 = pkin(3) * cos(t1026);
t1027 = qJ(2,1) + pkin(7);
t1079 = pkin(3) * cos(t1027);
t1078 = cos(pkin(7));
t1077 = 0.2e1 * pkin(2) * pkin(3);
t1022 = -pkin(6) - t1032;
t1016 = 0.1e1 / t1022;
t1036 = sin(qJ(1,3));
t1042 = cos(qJ(1,3));
t1035 = sin(qJ(2,3));
t1041 = cos(qJ(2,3));
t1028 = sin(pkin(7));
t987 = -m(3) * pkin(2) - mrSges(3,1) * t1078 + mrSges(3,2) * t1028 - mrSges(2,1);
t998 = t1028 * mrSges(3,1) + t1078 * mrSges(3,2) + mrSges(2,2);
t999 = mrSges(1,1) + (m(2) + m(3)) * pkin(1);
t1056 = t998 * t1035 + t987 * t1041 - t999;
t1066 = -m(2) * pkin(5) + mrSges(1,2) - mrSges(2,3) - mrSges(3,3);
t1059 = -t1032 * m(3) + t1066;
t1029 = legFrame(3,3);
t1010 = sin(t1029);
t1013 = cos(t1029);
t989 = -t1010 * g(1) + t1013 * g(2);
t992 = t1013 * g(1) + t1010 * g(2);
t968 = (-t1056 * t1036 + t1059 * t1042) * t992 + (t1059 * t1036 + t1056 * t1042) * t989;
t1076 = t1016 * t968;
t971 = t1036 * t992 - t1042 * t989;
t1075 = t1016 * t971;
t1019 = t1041 * pkin(2);
t995 = 0.1e1 / (t1019 + t1081);
t1074 = t1016 * t995;
t1023 = -pkin(6) - t1033;
t1017 = 0.1e1 / t1023;
t1038 = sin(qJ(1,2));
t1044 = cos(qJ(1,2));
t1037 = sin(qJ(2,2));
t1043 = cos(qJ(2,2));
t1055 = t998 * t1037 + t987 * t1043 - t999;
t1058 = -t1033 * m(3) + t1066;
t1030 = legFrame(2,3);
t1011 = sin(t1030);
t1014 = cos(t1030);
t990 = -t1011 * g(1) + t1014 * g(2);
t993 = t1014 * g(1) + t1011 * g(2);
t969 = (-t1055 * t1038 + t1058 * t1044) * t993 + (t1058 * t1038 + t1055 * t1044) * t990;
t1073 = t1017 * t969;
t972 = t1038 * t993 - t1044 * t990;
t1072 = t1017 * t972;
t1020 = t1043 * pkin(2);
t996 = 0.1e1 / (t1020 + t1080);
t1071 = t1017 * t996;
t1024 = -pkin(6) - t1034;
t1018 = 0.1e1 / t1024;
t1040 = sin(qJ(1,1));
t1046 = cos(qJ(1,1));
t1039 = sin(qJ(2,1));
t1045 = cos(qJ(2,1));
t1054 = t998 * t1039 + t987 * t1045 - t999;
t1057 = -t1034 * m(3) + t1066;
t1031 = legFrame(1,3);
t1012 = sin(t1031);
t1015 = cos(t1031);
t991 = -t1012 * g(1) + t1015 * g(2);
t994 = t1015 * g(1) + t1012 * g(2);
t970 = (-t1054 * t1040 + t1057 * t1046) * t994 + (t1057 * t1040 + t1054 * t1046) * t991;
t1070 = t1018 * t970;
t973 = t1040 * t994 - t1046 * t991;
t1069 = t1018 * t973;
t1021 = t1045 * pkin(2);
t997 = 0.1e1 / (t1021 + t1079);
t1068 = t1018 * t997;
t1067 = 0.2e1 * pkin(1);
t1065 = t1035 * pkin(2) + pkin(3) * sin(t1025);
t1064 = t1037 * pkin(2) + pkin(3) * sin(t1026);
t1063 = t1039 * pkin(2) + pkin(3) * sin(t1027);
t1062 = -t989 * t1036 - t992 * t1042;
t1061 = -t990 * t1038 - t993 * t1044;
t1060 = -t991 * t1040 - t994 * t1046;
t1053 = pkin(2) ^ 2;
t1052 = pkin(3) ^ 2;
t1051 = 0.2e1 * qJ(2,1);
t1050 = 0.2e1 * qJ(2,2);
t1049 = 0.2e1 * qJ(2,3);
t1002 = t1021 + pkin(1);
t1001 = t1020 + pkin(1);
t1000 = t1019 + pkin(1);
t988 = g(3) * t998;
t986 = g(3) * t987;
t985 = -t1012 * t1040 + t1015 * t1046;
t984 = t1012 * t1046 + t1015 * t1040;
t983 = -t1011 * t1038 + t1014 * t1044;
t982 = t1011 * t1044 + t1014 * t1038;
t981 = -t1010 * t1036 + t1013 * t1042;
t980 = t1010 * t1042 + t1013 * t1036;
t979 = t1002 * t1046 - t1040 * t1024;
t978 = t1001 * t1044 - t1038 * t1023;
t977 = t1000 * t1042 - t1036 * t1022;
t976 = t1040 * t1002 + t1046 * t1024;
t975 = t1038 * t1001 + t1044 * t1023;
t974 = t1036 * t1000 + t1042 * t1022;
t1 = [-t981 * t1076 - t983 * t1073 - t985 * t1070 - g(1) * m(4) + ((-t1012 * t976 + t979 * t1015 + t985 * t1079) * t1069 + (-t1011 * t975 + t978 * t1014 + t983 * t1080) * t1072 + (-t1010 * t974 + t977 * t1013 + t981 * t1081) * t1075) * m(3); -t980 * t1076 - t982 * t1073 - t984 * t1070 - g(2) * m(4) + ((t979 * t1012 + t976 * t1015 + t984 * t1079) * t1069 + (t978 * t1011 + t975 * t1014 + t982 * t1080) * t1072 + (t977 * t1010 + t974 * t1013 + t980 * t1081) * t1075) * m(3); -t1063 * t970 * t1068 + t997 * ((t1060 * t987 + t988) * t1039 - t1045 * (t1060 * t998 - t986)) - t1064 * t969 * t1071 + t996 * ((t1061 * t987 + t988) * t1037 - t1043 * (t1061 * t998 - t986)) - t1065 * t968 * t1074 + t995 * ((t1062 * t987 + t988) * t1035 - t1041 * (t1062 * t998 - t986)) - g(3) * m(4) + ((sin(t1051 + pkin(7)) * t1077 + t1052 * sin(0.2e1 * t1027) + t1053 * sin(t1051) + t1063 * t1067) * t973 * t1068 + (sin(t1050 + pkin(7)) * t1077 + t1052 * sin(0.2e1 * t1026) + t1053 * sin(t1050) + t1064 * t1067) * t972 * t1071 + (sin(t1049 + pkin(7)) * t1077 + t1052 * sin(0.2e1 * t1025) + t1053 * sin(t1049) + t1065 * t1067) * t971 * t1074) * m(3) / 0.2e1;];
taugX  = t1;
