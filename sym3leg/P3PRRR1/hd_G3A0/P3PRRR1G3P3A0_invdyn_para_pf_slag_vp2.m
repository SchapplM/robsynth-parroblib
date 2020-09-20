% Calculate vector of inverse dynamics forces for parallel robot
% P3PRRR1G3P3A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
% xDDP [3x1]
%   Generalized platform accelerations
% qJ [3x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% mrSges [4x3]
%   first moment of all robot links (mass times center of mass in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Ifges [4x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% tauX [3x1]
%   forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)
%   in platform coordinates xP

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:07
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3PRRR1G3P3A0_invdyn_para_pf_slag_vp2(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR1G3P3A0_invdyn_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRR1G3P3A0_invdyn_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3PRRR1G3P3A0_invdyn_para_pf_slag_vp2: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR1G3P3A0_invdyn_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRR1G3P3A0_invdyn_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3PRRR1G3P3A0_invdyn_para_pf_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRR1G3P3A0_invdyn_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3PRRR1G3P3A0_invdyn_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3PRRR1G3P3A0_invdyn_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR1G3P3A0_invdyn_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR1G3P3A0_invdyn_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:06:47
% EndTime: 2020-03-09 21:06:50
% DurationCPUTime: 2.46s
% Computational Cost: add. (15642->242), mult. (7887->437), div. (3366->5), fcn. (8976->24), ass. (0->187)
t133 = 0.1e1 / pkin(2);
t113 = pkin(7) + qJ(2,1);
t100 = cos(t113);
t103 = qJ(3,1) + t113;
t87 = sin(t103);
t90 = cos(t103);
t97 = sin(t113);
t214 = 0.1e1 / (t100 * t87 - t90 * t97);
t188 = t133 * t214;
t112 = pkin(7) + qJ(2,2);
t102 = qJ(3,2) + t112;
t86 = sin(t102);
t89 = cos(t102);
t96 = sin(t112);
t99 = cos(t112);
t215 = 0.1e1 / (t86 * t99 - t89 * t96);
t218 = t133 * t215;
t111 = pkin(7) + qJ(2,3);
t101 = qJ(3,3) + t111;
t85 = sin(t101);
t88 = cos(t101);
t95 = sin(t111);
t98 = cos(t111);
t216 = 0.1e1 / (t85 * t98 - t88 * t95);
t217 = t133 * t216;
t117 = legFrame(1,2);
t109 = cos(t117);
t227 = t109 * t188;
t116 = legFrame(2,2);
t108 = cos(t116);
t226 = t108 * t218;
t115 = legFrame(3,2);
t107 = cos(t115);
t225 = t107 * t217;
t131 = 0.1e1 / pkin(3);
t69 = pkin(2) * t100 + pkin(3) * t90;
t191 = t131 * t69;
t132 = pkin(2) ^ 2;
t143 = m(3) * t132 + Ifges(2,3) + Ifges(3,3);
t126 = cos(qJ(3,1));
t213 = mrSges(3,1) * pkin(2);
t167 = t126 * t213;
t123 = sin(qJ(3,1));
t211 = mrSges(3,2) * pkin(2);
t94 = t123 * t211;
t63 = t143 - 0.2e1 * t94 + 0.2e1 * t167;
t75 = Ifges(3,3) - t94 + t167;
t224 = t75 * t191 - t63 * t90;
t68 = pkin(2) * t99 + pkin(3) * t89;
t192 = t131 * t68;
t125 = cos(qJ(3,2));
t168 = t125 * t213;
t122 = sin(qJ(3,2));
t93 = t122 * t211;
t62 = t143 - 0.2e1 * t93 + 0.2e1 * t168;
t74 = Ifges(3,3) - t93 + t168;
t223 = t74 * t192 - t62 * t89;
t67 = pkin(2) * t98 + pkin(3) * t88;
t193 = t131 * t67;
t124 = cos(qJ(3,3));
t169 = t124 * t213;
t121 = sin(qJ(3,3));
t92 = t121 * t211;
t61 = t143 - 0.2e1 * t92 + 0.2e1 * t169;
t73 = Ifges(3,3) - t92 + t169;
t222 = t73 * t193 - t61 * t88;
t221 = Ifges(3,3) * t191 - t75 * t90;
t220 = Ifges(3,3) * t192 - t74 * t89;
t219 = Ifges(3,3) * t193 - t73 * t88;
t212 = mrSges(3,1) * g(3);
t210 = mrSges(3,2) * g(3);
t209 = pkin(2) * (mrSges(3,1) * t121 + mrSges(3,2) * t124);
t208 = pkin(2) * (mrSges(3,1) * t122 + mrSges(3,2) * t125);
t207 = pkin(2) * (mrSges(3,1) * t123 + mrSges(3,2) * t126);
t206 = g(3) * mrSges(2,2);
t127 = xDP(3);
t129 = xDP(1);
t162 = t107 * t216 * t67;
t146 = t129 * t162;
t140 = t131 * t146;
t104 = sin(t115);
t128 = xDP(2);
t181 = t104 * t128;
t64 = pkin(2) * t95 + pkin(3) * t85;
t185 = t64 * t131;
t171 = t129 * t133;
t196 = t107 * t88;
t43 = t216 * t171 * t196;
t16 = t43 + (-t140 + ((-t85 + t185) * t127 + (-t88 + t193) * t181) * t216) * t133;
t170 = t131 * t133;
t19 = (-t146 + (t64 * t127 + t67 * t181) * t216) * t170;
t205 = t16 * t19;
t161 = t108 * t215 * t68;
t145 = t129 * t161;
t139 = t131 * t145;
t105 = sin(t116);
t179 = t105 * t128;
t65 = pkin(2) * t96 + pkin(3) * t86;
t184 = t65 * t131;
t195 = t108 * t89;
t44 = t215 * t171 * t195;
t17 = t44 + (-t139 + ((-t86 + t184) * t127 + (-t89 + t192) * t179) * t215) * t133;
t20 = (-t145 + (t127 * t65 + t68 * t179) * t215) * t170;
t204 = t17 * t20;
t160 = t109 * t214 * t69;
t144 = t129 * t160;
t138 = t131 * t144;
t106 = sin(t117);
t177 = t106 * t128;
t66 = pkin(2) * t97 + pkin(3) * t87;
t183 = t66 * t131;
t194 = t109 * t90;
t45 = t214 * t171 * t194;
t18 = t45 + (-t138 + ((-t87 + t183) * t127 + (-t90 + t191) * t177) * t214) * t133;
t21 = (-t144 + (t66 * t127 + t69 * t177) * t214) * t170;
t203 = t18 * t21;
t119 = xDDP(2);
t182 = t104 * t119;
t180 = t105 * t119;
t178 = t106 * t119;
t120 = xDDP(1);
t174 = t109 * t120;
t172 = t120 * t131;
t114 = m(1) + m(2) + m(3);
t166 = (t104 * t107 + t105 * t108 + t106 * t109) * t114;
t156 = t104 * t217;
t155 = t105 * t218;
t154 = t106 * t188;
t153 = t120 * t196;
t152 = t120 * t195;
t151 = t90 * t174;
t150 = 0.2e1 * pkin(2) * pkin(3);
t79 = g(1) * t107 - g(2) * t104;
t149 = (mrSges(3,2) * t79 + t212) * t88 + (mrSges(3,1) * t79 - t210) * t85;
t80 = g(1) * t108 - g(2) * t105;
t148 = (mrSges(3,2) * t80 + t212) * t89 + (mrSges(3,1) * t80 - t210) * t86;
t81 = g(1) * t109 - g(2) * t106;
t147 = (mrSges(3,2) * t81 + t212) * t90 + (mrSges(3,1) * t81 - t210) * t87;
t142 = t85 * t95 + t88 * t98;
t141 = t86 * t96 + t89 * t99;
t137 = t100 * t90 + t87 * t97;
t136 = t142 * pkin(2);
t135 = t141 * pkin(2);
t134 = t137 * pkin(2);
t130 = pkin(3) ^ 2;
t118 = xDDP(3);
t110 = m(3) * pkin(2) + mrSges(2,1);
t91 = g(3) * t110;
t78 = g(1) * t106 + g(2) * t109;
t77 = g(1) * t105 + g(2) * t108;
t76 = g(1) * t104 + g(2) * t107;
t42 = (Ifges(3,3) * t183 - t75 * t87) * t188;
t41 = (Ifges(3,3) * t184 - t74 * t86) * t218;
t40 = (Ifges(3,3) * t185 - t73 * t85) * t217;
t39 = t221 * t227;
t38 = t220 * t226;
t37 = t219 * t225;
t36 = t221 * t154;
t35 = t220 * t155;
t34 = t219 * t156;
t33 = (t75 * t183 - t63 * t87) * t188;
t32 = (t74 * t184 - t62 * t86) * t218;
t31 = (t73 * t185 - t61 * t85) * t217;
t29 = t223 * t226;
t28 = t222 * t225;
t27 = t224 * t154;
t26 = t223 * t155;
t25 = t222 * t156;
t24 = t45 + (-t87 * t127 - t90 * t177) * t188;
t23 = t44 + (-t86 * t127 - t89 * t179) * t218;
t22 = t43 + (-t85 * t127 - t88 * t181) * t217;
t15 = t45 + (-t138 / 0.2e1 + ((-t87 + t183 / 0.2e1) * t127 + (-t90 + t191 / 0.2e1) * t177) * t214) * t133;
t14 = t44 + (-t139 / 0.2e1 + ((-t86 + t184 / 0.2e1) * t127 + (-t89 + t192 / 0.2e1) * t179) * t215) * t133;
t13 = t43 + (-t140 / 0.2e1 + ((-t85 + t185 / 0.2e1) * t127 + (-t88 + t193 / 0.2e1) * t181) * t216) * t133;
t12 = (pkin(3) * t203 + (pkin(3) * t18 + t24 * t134) * t24) * t188;
t11 = (pkin(3) * t204 + (t17 * pkin(3) + t23 * t135) * t23) * t218;
t10 = (pkin(3) * t205 + (t16 * pkin(3) + t22 * t136) * t22) * t217;
t9 = (-(-t137 * t15 * t150 - t130 * t18 - t132 * t24) * t131 * t24 + (pkin(3) + t134) * t203) * t188;
t8 = ((-t141 * t14 * t150 - t130 * t17 - t132 * t23) * t131 * t23 - (pkin(3) + t135) * t204) * t218;
t7 = ((-t142 * t13 * t150 - t130 * t16 - t132 * t22) * t131 * t22 - (pkin(3) + t136) * t205) * t217;
t6 = t24 ^ 2 * t207 - Ifges(3,3) * t9 + t12 * t75 + t147;
t5 = t23 ^ 2 * t208 + Ifges(3,3) * t8 + t11 * t74 + t148;
t4 = t22 ^ 2 * t209 + Ifges(3,3) * t7 + t10 * t73 + t149;
t3 = t63 * t12 - t75 * t9 - 0.2e1 * t21 * t15 * t207 + (mrSges(2,2) * t81 + t91) * t100 + t97 * (t110 * t81 - t206) + t147;
t2 = t62 * t11 + t74 * t8 - 0.2e1 * t20 * t14 * t208 + (mrSges(2,2) * t80 + t91) * t99 + t96 * (t110 * t80 - t206) + t148;
t1 = t61 * t10 + t73 * t7 - 0.2e1 * t19 * t13 * t209 + (mrSges(2,2) * t79 + t91) * t98 + t95 * (t110 * t79 - t206) + t149;
t30 = [(-g(1) + t120) * m(4) + t166 * t119 + (-t104 * t76 - t105 * t77 - t106 * t78 + (t104 ^ 2 + t105 ^ 2 + t106 ^ 2) * t120) * t114 + ((t3 * t194 - (-t118 * t87 + (t174 - t178) * t90) * t224 * t227) * t214 + (-t4 * t162 - t5 * t161 - t6 * t160 - (t118 * t66 + t69 * t178) * t214 * t39 + (t39 * t160 + t38 * t161 + t37 * t162) * t120) * t131 + ((-t38 * t184 + t29 * t86) * t118 + (-t120 * t29 + t2) * t195 + (-t38 * t192 + t29 * t89) * t180) * t215 + ((-t37 * t185 + t28 * t85) * t118 + (-t120 * t28 + t1) * t196 + (-t37 * t193 + t28 * t88) * t182) * t216) * t133; (-g(2) + t119) * m(4) + t166 * t120 + (-t107 * t76 - t108 * t77 - t109 * t78 + (t107 ^ 2 + t108 ^ 2 + t109 ^ 2) * t119) * t114 + ((-t36 * t160 - t35 * t161 - t34 * t162) * t172 + (t27 * t151 + (t36 * t183 - t27 * t87) * t118 + ((t36 * t191 - t27 * t90) * t119 - t90 * t3 + t6 * t191) * t106) * t214 + (t26 * t152 + (t35 * t184 - t26 * t86) * t118 + ((t35 * t192 - t26 * t89) * t119 - t89 * t2 + t5 * t192) * t105) * t215 + (t25 * t153 + (t34 * t185 - t25 * t85) * t118 + ((t34 * t193 - t25 * t88) * t119 - t88 * t1 + t4 * t193) * t104) * t216) * t133; (-g(3) + t118) * m(4) + ((-t42 * t160 - t41 * t161 - t40 * t162) * t172 + (t33 * t151 + (t42 * t183 - t33 * t87) * t118 - t87 * t3 + t6 * t183 + (t42 * t191 - t33 * t90) * t178) * t214 + (t32 * t152 + (t41 * t184 - t32 * t86) * t118 - t86 * t2 + t5 * t184 + (t41 * t192 - t32 * t89) * t180) * t215 + (t31 * t153 + (t40 * t185 - t31 * t85) * t118 - t85 * t1 + t4 * t185 + (t40 * t193 - t31 * t88) * t182) * t216) * t133;];
tauX  = t30;
