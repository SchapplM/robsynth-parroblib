% Calculate vector of inverse dynamics forces for parallel robot
% P3PRRR1G2P2A0
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
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Icges [4x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
%
% Output:
% tauX [3x1]
%   forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)
%   in platform coordinates xP

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:18
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3PRRR1G2P2A0_invdyn_para_pf_slag_vp1(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR1G2P2A0_invdyn_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRR1G2P2A0_invdyn_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3PRRR1G2P2A0_invdyn_para_pf_slag_vp1: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR1G2P2A0_invdyn_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRR1G2P2A0_invdyn_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3PRRR1G2P2A0_invdyn_para_pf_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRR1G2P2A0_invdyn_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3PRRR1G2P2A0_invdyn_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3PRRR1G2P2A0_invdyn_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR1G2P2A0_invdyn_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR1G2P2A0_invdyn_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:18:15
% EndTime: 2020-03-09 21:18:17
% DurationCPUTime: 2.18s
% Computational Cost: add. (17271->251), mult. (9210->479), div. (3366->5), fcn. (9498->18), ass. (0->200)
t128 = pkin(7) + qJ(2,3);
t117 = qJ(3,3) + t128;
t104 = sin(t117);
t107 = cos(t117);
t111 = sin(t128);
t114 = cos(t128);
t238 = 0.1e1 / (-t104 * t114 + t107 * t111);
t129 = pkin(7) + qJ(2,2);
t118 = qJ(3,2) + t129;
t105 = sin(t118);
t108 = cos(t118);
t112 = sin(t129);
t115 = cos(t129);
t237 = 0.1e1 / (-t105 * t115 + t108 * t112);
t130 = pkin(7) + qJ(2,1);
t119 = qJ(3,1) + t130;
t106 = sin(t119);
t109 = cos(t119);
t113 = sin(t130);
t116 = cos(t130);
t236 = 0.1e1 / (-t106 * t116 + t109 * t113);
t149 = 0.1e1 / pkin(3);
t151 = 0.1e1 / pkin(2);
t197 = t149 * t151;
t171 = t238 * t197;
t133 = legFrame(3,2);
t123 = cos(t133);
t145 = xDP(1);
t205 = t123 * t145;
t120 = sin(t133);
t144 = xDP(2);
t208 = t120 * t144;
t143 = xDP(3);
t76 = pkin(2) * t114 + pkin(3) * t107;
t211 = t76 * t143;
t73 = pkin(2) * t111 + pkin(3) * t104;
t34 = (t211 + (t205 - t208) * t73) * t171;
t214 = t151 * t238;
t174 = t120 * t214;
t229 = t104 * t238;
t179 = t123 * t229;
t198 = t145 * t151;
t199 = t143 * t151;
t226 = t107 * t238;
t40 = t104 * t144 * t174 - t179 * t198 - t199 * t226;
t16 = t34 + t40;
t235 = t16 * t34;
t170 = t237 * t197;
t134 = legFrame(2,2);
t124 = cos(t134);
t203 = t124 * t145;
t121 = sin(t134);
t207 = t121 * t144;
t77 = pkin(2) * t115 + pkin(3) * t108;
t210 = t77 * t143;
t74 = pkin(2) * t112 + pkin(3) * t105;
t35 = (t210 + (t203 - t207) * t74) * t170;
t213 = t151 * t237;
t173 = t121 * t213;
t228 = t105 * t237;
t177 = t124 * t228;
t225 = t108 * t237;
t41 = t105 * t144 * t173 - t177 * t198 - t199 * t225;
t17 = t35 + t41;
t234 = t17 * t35;
t169 = t236 * t197;
t135 = legFrame(1,2);
t125 = cos(t135);
t201 = t125 * t145;
t122 = sin(t135);
t206 = t122 * t144;
t78 = pkin(2) * t116 + pkin(3) * t109;
t209 = t78 * t143;
t75 = pkin(2) * t113 + pkin(3) * t106;
t36 = (t209 + (t201 - t206) * t75) * t169;
t212 = t151 * t236;
t172 = t122 * t212;
t227 = t106 * t236;
t175 = t125 * t227;
t224 = t109 * t236;
t42 = t106 * t144 * t172 - t175 * t198 - t199 * t224;
t18 = t36 + t42;
t233 = t18 * t36;
t232 = t238 * t76;
t231 = t237 * t77;
t230 = t236 * t78;
t223 = t120 * t238;
t222 = t121 * t237;
t221 = t122 * t236;
t220 = t149 * t238;
t219 = t149 * t237;
t218 = t149 * t73;
t217 = t149 * t74;
t216 = t149 * t75;
t126 = rSges(3,1) ^ 2 + rSges(3,2) ^ 2;
t94 = m(3) * t126 + Icges(3,3);
t215 = t149 * t94;
t204 = t123 * t151;
t202 = t124 * t151;
t200 = t125 * t151;
t131 = m(1) + m(2) + m(3);
t196 = (t120 * t123 + t121 * t124 + t122 * t125) * t131;
t195 = t73 * t223;
t194 = t74 * t222;
t193 = t75 * t221;
t192 = t123 * t238 * t73;
t191 = t124 * t237 * t74;
t190 = t125 * t236 * t75;
t189 = t238 * t218;
t188 = t76 * t220;
t187 = t237 * t217;
t186 = t77 * t219;
t185 = t236 * t216;
t184 = t149 * t230;
t183 = t73 * t215;
t182 = t74 * t215;
t181 = t75 * t215;
t180 = t104 * t223;
t178 = t105 * t222;
t176 = t106 * t221;
t168 = 0.2e1 * pkin(2) * pkin(3);
t167 = (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) + Icges(2,3) + Icges(3,3);
t139 = rSges(3,2) * g(3);
t141 = rSges(3,1) * g(3);
t91 = g(1) * t123 - g(2) * t120;
t166 = (rSges(3,2) * t91 + t141) * t104 - t107 * (rSges(3,1) * t91 - t139);
t79 = -rSges(3,1) * t114 - t111 * rSges(3,2);
t80 = rSges(3,1) * t111 - rSges(3,2) * t114;
t165 = t104 * t79 + t107 * t80;
t92 = g(1) * t124 - g(2) * t121;
t164 = (rSges(3,2) * t92 + t141) * t105 - t108 * (rSges(3,1) * t92 - t139);
t81 = -rSges(3,1) * t115 - t112 * rSges(3,2);
t82 = rSges(3,1) * t112 - rSges(3,2) * t115;
t163 = t105 * t81 + t108 * t82;
t93 = g(1) * t125 - g(2) * t122;
t162 = (rSges(3,2) * t93 + t141) * t106 - t109 * (rSges(3,1) * t93 - t139);
t83 = -rSges(3,1) * t116 - t113 * rSges(3,2);
t84 = rSges(3,1) * t113 - rSges(3,2) * t116;
t161 = t106 * t83 + t109 * t84;
t160 = t104 * t111 + t107 * t114;
t159 = t105 * t112 + t108 * t115;
t158 = t106 * t113 + t109 * t116;
t157 = (t104 * t80 - t107 * t79) * pkin(2);
t156 = (t105 * t82 - t108 * t81) * pkin(2);
t155 = (t106 * t84 - t109 * t83) * pkin(2);
t154 = t160 * pkin(2);
t153 = t159 * pkin(2);
t152 = t158 * pkin(2);
t150 = pkin(2) ^ 2;
t148 = pkin(3) ^ 2;
t142 = rSges(2,1) * g(3);
t140 = rSges(2,2) * g(3);
t138 = xDDP(1);
t137 = xDDP(2);
t136 = xDDP(3);
t132 = g(3) * pkin(2) * m(3);
t110 = t150 + t126;
t90 = g(1) * t122 + g(2) * t125;
t89 = g(1) * t121 + g(2) * t124;
t88 = g(1) * t120 + g(2) * t123;
t48 = Icges(3,3) + (t126 + t155) * m(3);
t47 = Icges(3,3) + (t126 + t156) * m(3);
t46 = Icges(3,3) + (t126 + t157) * m(3);
t45 = (t110 + 0.2e1 * t155) * m(3) + t167;
t44 = (t110 + 0.2e1 * t156) * m(3) + t167;
t43 = (t110 + 0.2e1 * t157) * m(3) + t167;
t39 = (t94 * t184 - t48 * t224) * t151;
t38 = (t94 * t186 - t47 * t225) * t151;
t37 = (t94 * t188 - t46 * t226) * t151;
t33 = (t106 * t48 - t181) * t172;
t32 = (t105 * t47 - t182) * t173;
t31 = (t104 * t46 - t183) * t174;
t30 = (t181 * t236 - t48 * t227) * t200;
t29 = (t182 * t237 - t47 * t228) * t202;
t28 = (t183 * t238 - t46 * t229) * t204;
t27 = (t48 * t184 - t45 * t224) * t151;
t26 = (t47 * t186 - t44 * t225) * t151;
t25 = (t46 * t188 - t43 * t226) * t151;
t24 = (t106 * t45 - t48 * t216) * t172;
t23 = (t105 * t44 - t47 * t217) * t173;
t22 = (t104 * t43 - t46 * t218) * t174;
t21 = (t48 * t185 - t45 * t227) * t200;
t20 = (t47 * t187 - t44 * t228) * t202;
t19 = (t46 * t189 - t43 * t229) * t204;
t15 = (t209 / 0.2e1 + (t201 / 0.2e1 - t206 / 0.2e1) * t75) * t169 + t42;
t14 = (t210 / 0.2e1 + (t203 / 0.2e1 - t207 / 0.2e1) * t74) * t170 + t41;
t13 = (t211 / 0.2e1 + (t205 / 0.2e1 - t208 / 0.2e1) * t73) * t171 + t40;
t12 = (pkin(3) * t233 + (t18 * pkin(3) + t42 * t152) * t42) * t212;
t11 = (pkin(3) * t234 + (t17 * pkin(3) + t41 * t153) * t41) * t213;
t10 = (pkin(3) * t235 + (t16 * pkin(3) + t40 * t154) * t40) * t214;
t9 = ((-t158 * t15 * t168 - t148 * t18 - t150 * t42) * t149 * t42 - (pkin(3) + t152) * t233) * t212;
t8 = ((-t159 * t14 * t168 - t148 * t17 - t150 * t41) * t41 * t219 - (pkin(3) + t153) * t237 * t234) * t151;
t7 = ((-t160 * t13 * t168 - t148 * t16 - t150 * t40) * t40 * t220 - (pkin(3) + t154) * t238 * t235) * t151;
t6 = -t48 * t12 - t94 * t9 + (-t161 * pkin(2) * t42 ^ 2 + t162) * m(3);
t5 = -t47 * t11 - t94 * t8 + (-t163 * pkin(2) * t41 ^ 2 + t164) * m(3);
t4 = -t46 * t10 - t94 * t7 + (-t165 * pkin(2) * t40 ^ 2 + t166) * m(3);
t3 = -t45 * t12 - t48 * t9 + m(2) * (-rSges(2,1) * t93 + t140) * t116 + t113 * (t132 + m(2) * (rSges(2,2) * t93 + t142)) + ((0.2e1 * t15 * t161 * t36 - t93 * t116) * pkin(2) + t162) * m(3);
t2 = -t44 * t11 - t47 * t8 + m(2) * (-rSges(2,1) * t92 + t140) * t115 + t112 * (t132 + m(2) * (rSges(2,2) * t92 + t142)) + ((0.2e1 * t14 * t163 * t35 - t92 * t115) * pkin(2) + t164) * m(3);
t1 = -t43 * t10 - t46 * t7 + m(2) * (-rSges(2,1) * t91 + t140) * t114 + t111 * (t132 + m(2) * (rSges(2,2) * t91 + t142)) + ((0.2e1 * t13 * t165 * t34 - t91 * t114) * pkin(2) + t166) * m(3);
t49 = [(-g(1) + t138) * m(4) + t196 * t137 + (-t120 * t88 - t121 * t89 - t122 * t90 + (t120 ^ 2 + t121 ^ 2 + t122 ^ 2) * t138) * t131 + ((t19 * t180 + t20 * t178 + t21 * t176 + (-t30 * t193 - t29 * t194 - t28 * t195) * t149) * t137 + (-t19 * t226 - t20 * t225 - t21 * t224 + (t30 * t230 + t29 * t231 + t28 * t232) * t149) * t136 + ((t30 * t185 - t21 * t227) * t138 - t3 * t227 + t6 * t185) * t125 + ((t29 * t187 - t20 * t228) * t138 - t2 * t228 + t5 * t187) * t124 + ((t28 * t189 - t19 * t229) * t138 - t1 * t229 + t4 * t189) * t123) * t151; (-g(2) + t137) * m(4) + t196 * t138 + (-t123 * t88 - t124 * t89 - t125 * t90 + (t123 ^ 2 + t124 ^ 2 + t125 ^ 2) * t137) * t131 + ((-t22 * t179 - t23 * t177 - t24 * t175 + (t33 * t190 + t32 * t191 + t31 * t192) * t149) * t138 + (-t22 * t226 - t23 * t225 - t24 * t224 + (t33 * t230 + t32 * t231 + t31 * t232) * t149) * t136 + ((t106 * t24 - t33 * t216) * t137 + t106 * t3 - t6 * t216) * t221 + ((t105 * t23 - t32 * t217) * t137 + t105 * t2 - t5 * t217) * t222 + ((t104 * t22 - t31 * t218) * t137 + t104 * t1 - t4 * t218) * t223) * t151; (-g(3) + t136) * m(4) + (-t1 * t226 - t2 * t225 - t3 * t224 + (t6 * t230 + t5 * t231 + t4 * t232) * t149 + (-t25 * t179 - t26 * t177 - t27 * t175 + (t39 * t190 + t38 * t191 + t37 * t192) * t149) * t138 + (t25 * t180 + t26 * t178 + t27 * t176 + (-t39 * t193 - t38 * t194 - t37 * t195) * t149) * t137 + (-t25 * t226 - t26 * t225 - t27 * t224 + (t39 * t230 + t38 * t231 + t37 * t232) * t149) * t136) * t151;];
tauX  = t49;
