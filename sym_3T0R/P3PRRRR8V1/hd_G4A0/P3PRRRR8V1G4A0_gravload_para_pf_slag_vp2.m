% Calculate Gravitation load for parallel robot
% P3PRRRR8V1G4A0
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4,alpha2,d2,d4,theta1]';
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
% Datum: 2020-08-06 17:27
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3PRRRR8V1G4A0_gravload_para_pf_slag_vp2(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V1G4A0_gravload_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V1G4A0_gravload_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3PRRRR8V1G4A0_gravload_para_pf_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR8V1G4A0_gravload_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRRR8V1G4A0_gravload_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3PRRRR8V1G4A0_gravload_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V1G4A0_gravload_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V1G4A0_gravload_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 17:25:21
% EndTime: 2020-08-06 17:25:23
% DurationCPUTime: 1.78s
% Computational Cost: add. (1623->206), mult. (3939->418), div. (54->7), fcn. (4263->34), ass. (0->165)
t1102 = sin(qJ(2,2));
t1108 = cos(qJ(2,2));
t1093 = legFrame(2,1);
t1067 = sin(t1093);
t1073 = cos(t1093);
t1097 = legFrame(2,2);
t1076 = sin(t1097);
t1079 = cos(t1097);
t1025 = t1076 * g(1) + (-g(2) * t1067 + g(3) * t1073) * t1079;
t1086 = sin(pkin(3));
t1088 = cos(pkin(3));
t1090 = legFrame(2,3);
t1064 = sin(t1090);
t1070 = cos(t1090);
t1149 = t1073 * t1076;
t1152 = t1067 * t1076;
t1163 = g(1) * t1079;
t1000 = t1070 * t1163 + (t1073 * t1064 + t1070 * t1152) * g(2) + (t1067 * t1064 - t1070 * t1149) * g(3);
t1085 = sin(pkin(6));
t1087 = cos(pkin(6));
t999 = -t1064 * t1163 + (-t1064 * t1152 + t1073 * t1070) * g(2) + (t1064 * t1149 + t1067 * t1070) * g(3);
t1128 = t1085 * t1000 - t1087 * t999;
t1126 = t1025 * t1086 - t1128 * t1088;
t1129 = t1000 * t1087 + t1085 * t999;
t1113 = t1126 * t1102 + t1129 * t1108;
t1100 = sin(qJ(2,3));
t1106 = cos(qJ(2,3));
t1092 = legFrame(3,1);
t1066 = sin(t1092);
t1072 = cos(t1092);
t1096 = legFrame(3,2);
t1075 = sin(t1096);
t1078 = cos(t1096);
t1024 = t1075 * g(1) + (-g(2) * t1066 + g(3) * t1072) * t1078;
t1089 = legFrame(3,3);
t1063 = sin(t1089);
t1069 = cos(t1089);
t1150 = t1072 * t1075;
t1153 = t1066 * t1075;
t1164 = g(1) * t1078;
t997 = -t1063 * t1164 + (-t1063 * t1153 + t1072 * t1069) * g(2) + (t1063 * t1150 + t1066 * t1069) * g(3);
t998 = t1069 * t1164 + (t1072 * t1063 + t1069 * t1153) * g(2) + (t1066 * t1063 - t1069 * t1150) * g(3);
t1131 = t1085 * t998 - t1087 * t997;
t1127 = t1024 * t1086 - t1131 * t1088;
t1130 = t1085 * t997 + t1087 * t998;
t1114 = t1127 * t1100 + t1130 * t1106;
t1105 = cos(qJ(3,3));
t1167 = pkin(2) * t1105;
t1048 = -t1106 * pkin(5) + t1100 * t1167;
t1099 = sin(qJ(3,3));
t1170 = pkin(2) * t1099;
t1033 = t1048 * t1086 + t1088 * t1170;
t1173 = 0.1e1 / t1033;
t1181 = 0.1e1 / t1105 * t1173;
t1107 = cos(qJ(3,2));
t1166 = pkin(2) * t1107;
t1049 = -t1108 * pkin(5) + t1102 * t1166;
t1101 = sin(qJ(3,2));
t1169 = pkin(2) * t1101;
t1034 = t1049 * t1086 + t1088 * t1169;
t1172 = 0.1e1 / t1034;
t1180 = 0.1e1 / t1107 * t1172;
t1109 = cos(qJ(3,1));
t1104 = sin(qJ(2,1));
t1110 = cos(qJ(2,1));
t1165 = pkin(2) * t1109;
t1050 = -t1110 * pkin(5) + t1104 * t1165;
t1103 = sin(qJ(3,1));
t1168 = pkin(2) * t1103;
t1035 = t1050 * t1086 + t1088 * t1168;
t1171 = 0.1e1 / t1035;
t1179 = 0.1e1 / t1109 * t1171;
t1094 = legFrame(1,1);
t1068 = sin(t1094);
t1074 = cos(t1094);
t1098 = legFrame(1,2);
t1077 = sin(t1098);
t1080 = cos(t1098);
t1026 = t1077 * g(1) + (-g(2) * t1068 + g(3) * t1074) * t1080;
t1091 = legFrame(1,3);
t1065 = sin(t1091);
t1071 = cos(t1091);
t1148 = t1074 * t1077;
t1151 = t1068 * t1077;
t1162 = g(1) * t1080;
t1001 = -t1065 * t1162 + (-t1065 * t1151 + t1074 * t1071) * g(2) + (t1065 * t1148 + t1068 * t1071) * g(3);
t1002 = t1071 * t1162 + (t1074 * t1065 + t1071 * t1151) * g(2) + (t1068 * t1065 - t1071 * t1148) * g(3);
t993 = t1087 * t1001 - t1085 * t1002;
t1176 = -t1026 * t1088 + t993 * t1086;
t1175 = t1025 * t1088 + t1128 * t1086;
t1174 = t1024 * t1088 + t1131 * t1086;
t1124 = t1001 * t1085 + t1002 * t1087;
t1125 = t1026 * t1086 + t1088 * t993;
t1112 = t1125 * t1104 + t1124 * t1110;
t1161 = t1024 * t1173;
t1158 = t1025 * t1172;
t1155 = t1026 * t1171;
t1143 = t1088 * t1100;
t1042 = t1085 * t1106 + t1087 * t1143;
t1045 = -t1085 * t1143 + t1087 * t1106;
t1021 = t1042 * t1069 + t1063 * t1045;
t1147 = t1075 * t1021;
t1146 = t1086 * t1105;
t1145 = t1086 * t1107;
t1144 = t1086 * t1109;
t1142 = t1088 * t1102;
t1141 = t1088 * t1104;
t1140 = t1088 * t1106;
t1139 = t1088 * t1108;
t1138 = t1088 * t1110;
t1137 = ((-t1174 * mrSges(3,1) + t1114 * mrSges(3,2)) * t1105 + (t1114 * mrSges(3,1) + t1174 * mrSges(3,2)) * t1099) * t1181;
t1136 = ((-t1175 * mrSges(3,1) + t1113 * mrSges(3,2)) * t1107 + (t1113 * mrSges(3,1) + t1175 * mrSges(3,2)) * t1101) * t1180;
t1135 = ((t1176 * mrSges(3,1) + t1112 * mrSges(3,2)) * t1109 + (t1112 * mrSges(3,1) - t1176 * mrSges(3,2)) * t1103) * t1179;
t1095 = mrSges(2,2) - mrSges(3,3);
t1134 = (t1114 * t1095 + (t1130 * t1100 - t1127 * t1106) * (t1105 * mrSges(3,1) - t1099 * mrSges(3,2) + mrSges(2,1))) * t1181;
t1133 = (t1113 * t1095 + (t1129 * t1102 - t1126 * t1108) * (t1107 * mrSges(3,1) - t1101 * mrSges(3,2) + mrSges(2,1))) * t1180;
t1132 = (t1112 * t1095 + (t1124 * t1104 - t1125 * t1110) * (t1109 * mrSges(3,1) - t1103 * mrSges(3,2) + mrSges(2,1))) * t1179;
t1123 = t1063 * t1042 - t1045 * t1069;
t1043 = t1085 * t1108 + t1087 * t1142;
t1046 = -t1085 * t1142 + t1087 * t1108;
t1122 = t1064 * t1043 - t1046 * t1070;
t1044 = t1085 * t1110 + t1087 * t1141;
t1047 = -t1085 * t1141 + t1087 * t1110;
t1121 = t1065 * t1044 - t1047 * t1071;
t1120 = -t1048 * t1088 + t1086 * t1170;
t1119 = -t1049 * t1088 + t1086 * t1169;
t1118 = -t1050 * t1088 + t1086 * t1168;
t1111 = 0.1e1 / pkin(2);
t1081 = m(1) + m(2) + m(3);
t1053 = pkin(5) * t1104 + t1110 * t1165;
t1052 = pkin(5) * t1102 + t1108 * t1166;
t1051 = pkin(5) * t1100 + t1106 * t1167;
t1041 = -t1065 * t1085 + t1071 * t1087;
t1040 = -t1064 * t1085 + t1070 * t1087;
t1039 = -t1063 * t1085 + t1069 * t1087;
t1038 = t1065 * t1087 + t1071 * t1085;
t1037 = t1064 * t1087 + t1070 * t1085;
t1036 = t1063 * t1087 + t1069 * t1085;
t1023 = t1044 * t1071 + t1065 * t1047;
t1022 = t1043 * t1070 + t1064 * t1046;
t1020 = -t1085 * t1053 + t1118 * t1087;
t1019 = -t1085 * t1052 + t1119 * t1087;
t1018 = -t1085 * t1051 + t1120 * t1087;
t1017 = t1053 * t1087 + t1118 * t1085;
t1016 = t1052 * t1087 + t1119 * t1085;
t1015 = t1051 * t1087 + t1120 * t1085;
t1014 = t1038 * t1151 - t1041 * t1074;
t1013 = t1037 * t1152 - t1040 * t1073;
t1012 = t1036 * t1153 - t1039 * t1072;
t1011 = t1038 * t1074 + t1041 * t1151;
t1010 = t1037 * t1073 + t1040 * t1152;
t1009 = t1036 * t1072 + t1039 * t1153;
t1008 = -t1068 * t1038 + t1041 * t1148;
t1007 = t1038 * t1148 + t1068 * t1041;
t1006 = -t1067 * t1037 + t1040 * t1149;
t1005 = t1037 * t1149 + t1067 * t1040;
t1004 = -t1066 * t1036 + t1039 * t1150;
t1003 = t1036 * t1150 + t1066 * t1039;
t996 = -t1017 * t1065 + t1020 * t1071;
t995 = -t1016 * t1064 + t1019 * t1070;
t994 = -t1015 * t1063 + t1018 * t1069;
t990 = -t1080 * t1035 + (t1017 * t1071 + t1020 * t1065) * t1077;
t989 = -t1079 * t1034 + (t1016 * t1070 + t1019 * t1064) * t1076;
t988 = -t1078 * t1033 + (t1015 * t1069 + t1018 * t1063) * t1075;
t1 = [-(t1023 * t1103 + t1041 * t1144) * t1080 * t1132 - (t1022 * t1101 + t1040 * t1145) * t1079 * t1133 - (t1021 * t1099 + t1039 * t1146) * t1078 * t1134 - g(1) * m(4) + (-t1080 * ((-t1104 * t1038 + t1041 * t1138) * t1165 + pkin(5) * (t1110 * t1038 + t1041 * t1141)) * t1135 - t1079 * ((-t1102 * t1037 + t1040 * t1139) * t1166 + pkin(5) * (t1108 * t1037 + t1040 * t1142)) * t1136 - t1078 * ((-t1100 * t1036 + t1039 * t1140) * t1167 + pkin(5) * (t1106 * t1036 + t1039 * t1143)) * t1137) * t1111 + (-((t1118 * t1038 + t1041 * t1053) * t1080 + t1077 * t1035) * t1155 - ((t1119 * t1037 + t1040 * t1052) * t1079 + t1076 * t1034) * t1158 - ((t1120 * t1036 + t1039 * t1051) * t1078 + t1075 * t1033) * t1161) * t1081; ((-t1023 * t1151 - t1121 * t1074) * t1103 - t1011 * t1144) * t1132 + ((-t1022 * t1152 - t1122 * t1073) * t1101 - t1010 * t1145) * t1133 + ((-t1066 * t1147 - t1123 * t1072) * t1099 - t1009 * t1146) * t1134 - g(2) * m(4) + ((-(t1011 * t1138 - t1014 * t1104) * t1165 - pkin(5) * (t1011 * t1141 + t1110 * t1014)) * t1135 + (-(t1010 * t1139 - t1013 * t1102) * t1166 - pkin(5) * (t1010 * t1142 + t1108 * t1013)) * t1136 + (-(t1009 * t1140 - t1012 * t1100) * t1167 - pkin(5) * (t1009 * t1143 + t1106 * t1012)) * t1137) * t1111 + (-(t990 * t1068 - t1074 * t996) * t1155 - (t989 * t1067 - t1073 * t995) * t1158 - (t988 * t1066 - t1072 * t994) * t1161) * t1081; ((t1023 * t1148 - t1121 * t1068) * t1103 + t1008 * t1144) * t1132 + ((t1022 * t1149 - t1122 * t1067) * t1101 + t1006 * t1145) * t1133 + ((-t1123 * t1066 + t1072 * t1147) * t1099 + t1004 * t1146) * t1134 - g(3) * m(4) + (((-t1104 * t1007 + t1008 * t1138) * t1165 + pkin(5) * (t1110 * t1007 + t1008 * t1141)) * t1135 + ((-t1102 * t1005 + t1006 * t1139) * t1166 + pkin(5) * (t1108 * t1005 + t1006 * t1142)) * t1136 + ((-t1100 * t1003 + t1004 * t1140) * t1167 + pkin(5) * (t1106 * t1003 + t1004 * t1143)) * t1137) * t1111 + (-(-t1068 * t996 - t990 * t1074) * t1155 - (-t1067 * t995 - t989 * t1073) * t1158 - (-t1066 * t994 - t988 * t1072) * t1161) * t1081;];
taugX  = t1;
