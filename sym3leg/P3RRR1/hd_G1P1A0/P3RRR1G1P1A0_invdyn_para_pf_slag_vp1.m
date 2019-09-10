% Calculate vector of inverse dynamics forces for parallel robot
% P3RRR1G1P1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
% xDDP [3x1]
%   Generalized platform accelerations
% qJ [2x3]
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2,d3]';
% m [3x1]
%   mass of all robot links (leg links until cut joint, platform)
% rSges [3x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Icges [3x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
%
% Output:
% tauX [3x1]
%   forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)
%   in platform coordinates xP

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-03 15:38
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3RRR1G1P1A0_invdyn_para_pf_slag_vp1(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(2,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(5,1),zeros(2+1,1),zeros(2+1,3),zeros(2+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRR1G1P1A0_invdyn_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRR1G1P1A0_invdyn_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RRR1G1P1A0_invdyn_para_pf_slag_vp1: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [2 3]), ...
  'P3RRR1G1P1A0_invdyn_para_pf_slag_vp1: qJ has to be [2x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRR1G1P1A0_invdyn_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'P3RRR1G1P1A0_invdyn_para_pf_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [3 1]), ...
  'P3RRR1G1P1A0_invdyn_para_pf_slag_vp1: m has to be [3x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [3,3]), ...
  'P3RRR1G1P1A0_invdyn_para_pf_slag_vp1: rSges has to be [3x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [3 6]), ...
  'P3RRR1G1P1A0_invdyn_para_pf_slag_vp1: Icges has to be [3x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRR1G1P1A0_invdyn_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRR1G1P1A0_invdyn_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 15:38:12
% EndTime: 2019-05-03 15:38:14
% DurationCPUTime: 2.36s
% Computational Cost: add. (10265->250), mult. (11981->468), div. (2388->5), fcn. (12026->20), ass. (0->187)
t181 = 0.1e1 / pkin(2);
t183 = 0.1e1 / pkin(1);
t215 = t181 * t183;
t151 = qJ(1,3) + qJ(2,3);
t131 = sin(t151);
t134 = cos(t151);
t160 = sin(qJ(1,3));
t163 = cos(qJ(1,3));
t234 = 0.1e1 / (t131 * t163 - t134 * t160);
t207 = t234 * t215;
t154 = legFrame(3,3);
t137 = sin(t154);
t140 = cos(t154);
t94 = t131 * t140 + t134 * t137;
t67 = pkin(1) * (t137 * t163 + t140 * t160) + t94 * pkin(2);
t95 = -t131 * t137 + t134 * t140;
t70 = pkin(1) * (-t137 * t160 + t140 * t163) + t95 * pkin(2);
t169 = xP(3);
t148 = sin(t169);
t149 = cos(t169);
t172 = koppelP(3,2);
t175 = koppelP(3,1);
t114 = -t148 * t172 + t149 * t175;
t166 = xDP(3);
t167 = xDP(2);
t88 = t114 * t166 + t167;
t111 = t148 * t175 + t149 * t172;
t168 = xDP(1);
t91 = -t111 * t166 + t168;
t31 = (-t67 * t88 - t70 * t91) * t207;
t152 = qJ(1,2) + qJ(2,2);
t132 = sin(t152);
t135 = cos(t152);
t161 = sin(qJ(1,2));
t164 = cos(qJ(1,2));
t233 = 0.1e1 / (t132 * t164 - t135 * t161);
t205 = t233 * t215;
t155 = legFrame(2,3);
t138 = sin(t155);
t141 = cos(t155);
t96 = t132 * t141 + t135 * t138;
t68 = pkin(1) * (t138 * t164 + t141 * t161) + t96 * pkin(2);
t97 = -t132 * t138 + t135 * t141;
t71 = pkin(1) * (-t138 * t161 + t141 * t164) + t97 * pkin(2);
t173 = koppelP(2,2);
t176 = koppelP(2,1);
t115 = -t148 * t173 + t149 * t176;
t89 = t115 * t166 + t167;
t112 = t148 * t176 + t149 * t173;
t92 = -t112 * t166 + t168;
t32 = (-t68 * t89 - t71 * t92) * t205;
t153 = qJ(1,1) + qJ(2,1);
t133 = sin(t153);
t136 = cos(t153);
t162 = sin(qJ(1,1));
t165 = cos(qJ(1,1));
t232 = 0.1e1 / (t133 * t165 - t136 * t162);
t203 = t232 * t215;
t156 = legFrame(1,3);
t139 = sin(t156);
t142 = cos(t156);
t98 = t133 * t142 + t136 * t139;
t69 = pkin(1) * (t139 * t165 + t142 * t162) + t98 * pkin(2);
t99 = -t133 * t139 + t136 * t142;
t72 = pkin(1) * (-t139 * t162 + t142 * t165) + t99 * pkin(2);
t174 = koppelP(1,2);
t177 = koppelP(1,1);
t116 = -t148 * t174 + t149 * t177;
t90 = t116 * t166 + t167;
t113 = t148 * t177 + t149 * t174;
t93 = -t113 * t166 + t168;
t33 = (-t69 * t90 - t72 * t93) * t203;
t221 = t233 * t183;
t222 = t234 * t183;
t43 = (t88 * t94 + t91 * t95) * t222;
t237 = 0.2e1 * t43 + t31;
t44 = (t89 * t96 + t92 * t97) * t221;
t236 = 0.2e1 * t44 + t32;
t220 = t232 * t183;
t45 = (t90 * t98 + t93 * t99) * t220;
t235 = 0.2e1 * t45 + t33;
t184 = (t133 * t162 + t136 * t165) * pkin(1);
t185 = (t132 * t161 + t135 * t164) * pkin(1);
t186 = (t131 * t160 + t134 * t163) * pkin(1);
t22 = t43 + t31;
t231 = t22 * t31;
t23 = t44 + t32;
t230 = t23 * t32;
t24 = t45 + t33;
t229 = t24 * t33;
t228 = t234 * t94;
t227 = t234 * t95;
t226 = t233 * t96;
t225 = t233 * t97;
t224 = t232 * t98;
t223 = t232 * t99;
t219 = t234 * t181;
t218 = t233 * t181;
t217 = t232 * t181;
t146 = rSges(2,1) ^ 2 + rSges(2,2) ^ 2;
t129 = m(2) * t146 + Icges(2,3);
t216 = t129 * t181;
t214 = t67 * t219;
t213 = t70 * t219;
t212 = t68 * t218;
t211 = t71 * t218;
t210 = t69 * t217;
t209 = t72 * t217;
t208 = t234 * t216;
t206 = t233 * t216;
t204 = t232 * t216;
t202 = (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) + Icges(1,3) + Icges(2,3);
t117 = -g(1) * t137 + g(2) * t140;
t120 = g(1) * t140 + g(2) * t137;
t201 = t131 * (rSges(2,1) * t120 + rSges(2,2) * t117) - t134 * (rSges(2,1) * t117 - rSges(2,2) * t120);
t118 = -g(1) * t138 + g(2) * t141;
t121 = g(1) * t141 + g(2) * t138;
t200 = t132 * (rSges(2,1) * t121 + rSges(2,2) * t118) - t135 * (rSges(2,1) * t118 - rSges(2,2) * t121);
t119 = -g(1) * t139 + g(2) * t142;
t122 = g(1) * t142 + g(2) * t139;
t199 = t133 * (rSges(2,1) * t122 + rSges(2,2) * t119) - t136 * (rSges(2,1) * t119 - rSges(2,2) * t122);
t150 = t166 ^ 2;
t157 = xDDP(3);
t158 = xDDP(2);
t73 = -t111 * t150 + t114 * t157 + t158;
t159 = xDDP(1);
t76 = -t111 * t157 - t114 * t150 + t159;
t195 = t183 * (t67 * t73 + t70 * t76);
t74 = -t112 * t150 + t115 * t157 + t158;
t77 = -t112 * t157 - t115 * t150 + t159;
t194 = t183 * (t68 * t74 + t71 * t77);
t75 = -t113 * t150 + t116 * t157 + t158;
t78 = -t113 * t157 - t116 * t150 + t159;
t193 = t183 * (t69 * t75 + t72 * t78);
t192 = t183 * (t73 * t94 + t76 * t95);
t191 = t183 * (t74 * t96 + t77 * t97);
t190 = t183 * (t75 * t98 + t78 * t99);
t123 = -rSges(2,1) * t160 + rSges(2,2) * t163;
t126 = rSges(2,1) * t163 + rSges(2,2) * t160;
t189 = (-t123 * t131 + t126 * t134) * pkin(1);
t124 = -rSges(2,1) * t161 + rSges(2,2) * t164;
t127 = rSges(2,1) * t164 + rSges(2,2) * t161;
t188 = (-t124 * t132 + t127 * t135) * pkin(1);
t125 = -rSges(2,1) * t162 + rSges(2,2) * t165;
t128 = rSges(2,1) * t165 + rSges(2,2) * t162;
t187 = (-t125 * t133 + t128 * t136) * pkin(1);
t182 = pkin(1) ^ 2;
t180 = pkin(2) ^ 2;
t171 = rSges(3,1);
t170 = rSges(3,2);
t130 = t182 + t146;
t110 = -t148 * t170 + t149 * t171;
t109 = t148 * t171 + t149 * t170;
t81 = -t125 * t136 - t128 * t133;
t80 = -t124 * t135 - t127 * t132;
t79 = -t123 * t134 - t126 * t131;
t66 = Icges(2,3) + (t146 + t187) * m(2);
t65 = Icges(2,3) + (t146 + t188) * m(2);
t64 = Icges(2,3) + (t146 + t189) * m(2);
t63 = (t130 + 0.2e1 * t187) * m(2) + t202;
t62 = (t130 + 0.2e1 * t188) * m(2) + t202;
t61 = (t130 + 0.2e1 * t189) * m(2) + t202;
t48 = (-t113 * t99 + t116 * t98) * t220;
t47 = (-t112 * t97 + t115 * t96) * t221;
t46 = (-t111 * t95 + t114 * t94) * t222;
t42 = (-t113 * t72 + t116 * t69) * t203;
t41 = (-t112 * t71 + t115 * t68) * t205;
t40 = (-t111 * t70 + t114 * t67) * t207;
t18 = -t129 * t42 + t48 * t66;
t17 = -t129 * t41 + t47 * t65;
t16 = -t129 * t40 + t46 * t64;
t15 = -t42 * t66 + t48 * t63;
t14 = -t41 * t65 + t47 * t62;
t13 = -t40 * t64 + t46 * t61;
t12 = (-pkin(2) * t229 + (-pkin(2) * t24 - t45 * t184) * t45) * t220;
t11 = ((-pkin(2) * t23 - t44 * t185) * t44 - pkin(2) * t230) * t221;
t10 = ((-pkin(2) * t22 - t43 * t186) * t43 - pkin(2) * t231) * t222;
t9 = ((pkin(2) * t184 * t235 + t180 * t24 + t182 * t45) * t181 * t45 + (pkin(2) + t184) * t229) * t220;
t8 = ((pkin(2) * t185 * t236 + t180 * t23 + t182 * t44) * t181 * t44 + (pkin(2) + t185) * t230) * t221;
t7 = ((pkin(2) * t186 * t237 + t180 * t22 + t182 * t43) * t181 * t43 + (pkin(2) + t186) * t231) * t222;
t6 = -t66 * t12 - t129 * t9 + (-pkin(1) * t45 ^ 2 * t81 + t199) * m(2);
t5 = -t65 * t11 - t129 * t8 + (-pkin(1) * t44 ^ 2 * t80 + t200) * m(2);
t4 = -t64 * t10 - t129 * t7 + (-pkin(1) * t43 ^ 2 * t79 + t201) * m(2);
t3 = -t63 * t12 - t66 * t9 + ((-rSges(1,1) * t119 + rSges(1,2) * t122) * t165 + (rSges(1,1) * t122 + rSges(1,2) * t119) * t162) * m(1) + ((t33 * t81 * t235 - t119 * t165 + t122 * t162) * pkin(1) + t199) * m(2);
t2 = -t62 * t11 - t65 * t8 + ((-rSges(1,1) * t118 + rSges(1,2) * t121) * t164 + (rSges(1,1) * t121 + rSges(1,2) * t118) * t161) * m(1) + ((t32 * t80 * t236 - t118 * t164 + t121 * t161) * pkin(1) + t200) * m(2);
t1 = -t61 * t10 - t64 * t7 + ((-rSges(1,1) * t117 + rSges(1,2) * t120) * t163 + (rSges(1,1) * t120 + rSges(1,2) * t117) * t160) * m(1) + ((t31 * t79 * t237 - t117 * t163 + t120 * t160) * pkin(1) + t201) * m(2);
t19 = [(-t109 * t157 - t110 * t150 - g(1) + t159) * m(3) + ((t99 * t3 + (-t66 * t209 + t63 * t223) * t190) * t232 + (t97 * t2 + (-t65 * t211 + t62 * t225) * t191) * t233 + (t95 * t1 + (-t64 * t213 + t61 * t227) * t192) * t234 + (-(t6 * t72 + (-t72 * t204 + t66 * t223) * t193) * t232 - (t5 * t71 + (-t71 * t206 + t65 * t225) * t194) * t233 - (t4 * t70 + (-t70 * t208 + t64 * t227) * t195) * t234) * t181) * t183; (-t109 * t150 + t110 * t157 - g(2) + t158) * m(3) + ((t98 * t3 + (-t66 * t210 + t63 * t224) * t190) * t232 + (t96 * t2 + (-t65 * t212 + t62 * t226) * t191) * t233 + (t94 * t1 + (-t64 * t214 + t61 * t228) * t192) * t234 + (-(t6 * t69 + (-t69 * t204 + t66 * t224) * t193) * t232 - (t5 * t68 + (-t68 * t206 + t65 * t226) * t194) * t233 - (t4 * t67 + (-t67 * t208 + t64 * t228) * t195) * t234) * t181) * t183; Icges(3,3) * t157 + t46 * t1 + t47 * t2 + t48 * t3 - t40 * t4 - t41 * t5 - t42 * t6 + (-t109 * t159 + t110 * t158 + (t170 ^ 2 + t171 ^ 2) * t157 + (g(1) * t171 + g(2) * t170) * t148 + (g(1) * t170 - g(2) * t171) * t149) * m(3) + ((t15 * t223 - t18 * t209) * t78 + (t15 * t224 - t18 * t210) * t75 + (t14 * t225 - t17 * t211) * t77 + (t14 * t226 - t17 * t212) * t74 + (t13 * t227 - t16 * t213) * t76 + (t13 * t228 - t16 * t214) * t73) * t183;];
tauX  = t19;
