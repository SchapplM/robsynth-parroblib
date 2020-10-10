% Calculate Gravitation load for parallel robot
% P4PRRRR1G2A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [4x1]
%   Generalized platform coordinates
% qJ [3x4]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [4x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% koppelP [4x3]
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
% taugX [4x1]
%   forces required to compensate gravitation load
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 11:00
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P4PRRRR1G2A0_gravload_para_pf_slag_vp2(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,4),zeros(3,1),zeros(4,3),zeros(4,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR1G2A0_gravload_para_pf_slag_vp2: xP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR1G2A0_gravload_para_pf_slag_vp2: qJ has to be [3x4] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P4PRRRR1G2A0_gravload_para_pf_slag_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P4PRRRR1G2A0_gravload_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P4PRRRR1G2A0_gravload_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P4PRRRR1G2A0_gravload_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR1G2A0_gravload_para_pf_slag_vp2: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR1G2A0_gravload_para_pf_slag_vp2: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 10:58:50
% EndTime: 2020-08-07 10:58:51
% DurationCPUTime: 0.68s
% Computational Cost: add. (463->113), mult. (857->211), div. (128->13), fcn. (714->26), ass. (0->111)
t1116 = legFrame(1,2);
t1088 = sin(t1116);
t1092 = cos(t1116);
t1079 = t1088 * g(1) + t1092 * g(2);
t1112 = mrSges(2,2) - mrSges(3,3);
t1084 = t1112 * g(3);
t1122 = sin(qJ(2,1));
t1128 = cos(qJ(2,1));
t1129 = mrSges(2,1) * g(3);
t1121 = sin(qJ(3,1));
t1127 = cos(qJ(3,1));
t1189 = mrSges(3,1) * t1127 - t1121 * mrSges(3,2);
t1055 = t1084 * t1128 + t1122 * (t1189 * g(3) + t1129) + ((-mrSges(2,1) - t1189) * t1128 + t1122 * t1112) * t1079;
t1083 = t1092 * g(1) - t1088 * g(2);
t1106 = 0.1e1 / t1127;
t1154 = g(3) * t1128 + t1079 * t1122;
t1101 = 0.1e1 / t1122;
t1169 = t1101 * t1128;
t1142 = t1121 / t1127 ^ 2 * t1055 * t1169 - t1106 * ((mrSges(3,1) * t1083 + mrSges(3,2) * t1154) * t1127 + t1121 * (mrSges(3,1) * t1154 - mrSges(3,2) * t1083));
t1115 = legFrame(2,2);
t1087 = sin(t1115);
t1091 = cos(t1115);
t1078 = t1087 * g(1) + t1091 * g(2);
t1120 = sin(qJ(2,2));
t1126 = cos(qJ(2,2));
t1119 = sin(qJ(3,2));
t1125 = cos(qJ(3,2));
t1188 = mrSges(3,1) * t1125 - t1119 * mrSges(3,2);
t1054 = t1084 * t1126 + t1120 * (t1188 * g(3) + t1129) + ((-mrSges(2,1) - t1188) * t1126 + t1120 * t1112) * t1078;
t1082 = t1091 * g(1) - t1087 * g(2);
t1104 = 0.1e1 / t1125;
t1155 = g(3) * t1126 + t1078 * t1120;
t1100 = 0.1e1 / t1120;
t1171 = t1100 * t1126;
t1143 = t1119 / t1125 ^ 2 * t1054 * t1171 - t1104 * ((mrSges(3,1) * t1082 + mrSges(3,2) * t1155) * t1125 + t1119 * (mrSges(3,1) * t1155 - mrSges(3,2) * t1082));
t1114 = legFrame(3,2);
t1086 = sin(t1114);
t1090 = cos(t1114);
t1077 = t1086 * g(1) + t1090 * g(2);
t1118 = sin(qJ(2,3));
t1124 = cos(qJ(2,3));
t1117 = sin(qJ(3,3));
t1123 = cos(qJ(3,3));
t1187 = mrSges(3,1) * t1123 - t1117 * mrSges(3,2);
t1053 = t1084 * t1124 + t1118 * (t1187 * g(3) + t1129) + ((-mrSges(2,1) - t1187) * t1124 + t1118 * t1112) * t1077;
t1081 = t1090 * g(1) - t1086 * g(2);
t1102 = 0.1e1 / t1123;
t1156 = g(3) * t1124 + t1077 * t1118;
t1099 = 0.1e1 / t1118;
t1173 = t1099 * t1124;
t1144 = t1117 / t1123 ^ 2 * t1053 * t1173 - t1102 * ((mrSges(3,1) * t1081 + mrSges(3,2) * t1156) * t1123 + t1117 * (mrSges(3,1) * t1156 - mrSges(3,2) * t1081));
t1113 = legFrame(4,2);
t1085 = sin(t1113);
t1089 = cos(t1113);
t1076 = t1085 * g(1) + t1089 * g(2);
t1109 = sin(qJ(2,4));
t1111 = cos(qJ(2,4));
t1108 = sin(qJ(3,4));
t1110 = cos(qJ(3,4));
t1186 = mrSges(3,1) * t1110 - t1108 * mrSges(3,2);
t1052 = t1084 * t1111 + t1109 * (t1186 * g(3) + t1129) + ((-mrSges(2,1) - t1186) * t1111 + t1109 * t1112) * t1076;
t1080 = t1089 * g(1) - t1085 * g(2);
t1096 = 0.1e1 / t1110;
t1157 = g(3) * t1111 + t1076 * t1109;
t1095 = 0.1e1 / t1109;
t1176 = t1095 * t1111;
t1145 = t1108 / t1110 ^ 2 * t1052 * t1176 - t1096 * ((mrSges(3,1) * t1080 + mrSges(3,2) * t1157) * t1110 + t1108 * (mrSges(3,1) * t1157 - mrSges(3,2) * t1080));
t1177 = t1095 * t1096;
t1174 = t1099 * t1102;
t1172 = t1100 * t1104;
t1170 = t1101 * t1106;
t1165 = t1109 * t1110;
t1164 = t1118 * t1123;
t1163 = t1120 * t1125;
t1162 = t1122 * t1127;
t1161 = t1076 * t1177;
t1160 = t1077 * t1174;
t1159 = t1078 * t1172;
t1158 = t1079 * t1170;
t1141 = 0.1e1 / pkin(2);
t1140 = koppelP(1,1);
t1139 = koppelP(2,1);
t1138 = koppelP(3,1);
t1137 = koppelP(4,1);
t1136 = koppelP(1,2);
t1135 = koppelP(2,2);
t1134 = koppelP(3,2);
t1133 = koppelP(4,2);
t1132 = mrSges(4,1);
t1131 = mrSges(4,2);
t1130 = xP(4);
t1098 = m(1) + m(2) + m(3);
t1094 = cos(t1130);
t1093 = sin(t1130);
t1075 = -t1093 * t1136 + t1094 * t1140;
t1074 = -t1093 * t1135 + t1094 * t1139;
t1073 = -t1093 * t1134 + t1094 * t1138;
t1072 = -t1093 * t1133 + t1094 * t1137;
t1071 = -t1093 * t1140 - t1094 * t1136;
t1070 = -t1093 * t1139 - t1094 * t1135;
t1069 = -t1093 * t1138 - t1094 * t1134;
t1068 = -t1093 * t1137 - t1094 * t1133;
t1067 = t1088 * t1121 + t1092 * t1162;
t1066 = t1087 * t1119 + t1091 * t1163;
t1065 = t1086 * t1117 + t1090 * t1164;
t1064 = t1088 * t1162 - t1092 * t1121;
t1063 = t1087 * t1163 - t1091 * t1119;
t1062 = t1086 * t1164 - t1090 * t1117;
t1061 = t1085 * t1108 + t1089 * t1165;
t1060 = t1085 * t1165 - t1089 * t1108;
t1 = [-g(1) * m(4) + (-t1060 * t1161 - t1062 * t1160 - t1063 * t1159 - t1064 * t1158) * t1098 + (t1089 * t1145 + t1090 * t1144 + t1091 * t1143 + t1092 * t1142) * t1141; -g(2) * m(4) + (-t1061 * t1161 - t1065 * t1160 - t1066 * t1159 - t1067 * t1158) * t1098 + (-t1085 * t1145 - t1086 * t1144 - t1087 * t1143 - t1088 * t1142) * t1141; -g(3) * m(4) + (-t1052 * t1177 - t1053 * t1174 - t1054 * t1172 - t1055 * t1170) * t1141 + (-t1076 * t1176 - t1077 * t1173 - t1078 * t1171 - t1079 * t1169) * t1098; -(-g(1) * t1132 - g(2) * t1131) * t1093 + t1094 * (g(1) * t1131 - g(2) * t1132) + (-(t1064 * t1071 + t1067 * t1075) * t1158 - (t1063 * t1070 + t1066 * t1074) * t1159 - (t1062 * t1069 + t1065 * t1073) * t1160 - (t1060 * t1068 + t1061 * t1072) * t1161) * t1098 + (t1145 * (t1068 * t1089 - t1072 * t1085) + t1144 * (t1069 * t1090 - t1073 * t1086) + t1143 * (t1070 * t1091 - t1074 * t1087) + t1142 * (t1071 * t1092 - t1075 * t1088)) * t1141;];
taugX  = t1;
