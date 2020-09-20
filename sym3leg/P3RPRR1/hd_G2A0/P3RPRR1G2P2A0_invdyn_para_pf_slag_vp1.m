% Calculate vector of inverse dynamics forces for parallel robot
% P3RPRR1G2P2A0
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
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2020-03-09 21:25
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3RPRR1G2P2A0_invdyn_para_pf_slag_vp1(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRR1G2P2A0_invdyn_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRR1G2P2A0_invdyn_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RPRR1G2P2A0_invdyn_para_pf_slag_vp1: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRR1G2P2A0_invdyn_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRR1G2P2A0_invdyn_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRR1G2P2A0_invdyn_para_pf_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRR1G2P2A0_invdyn_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RPRR1G2P2A0_invdyn_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RPRR1G2P2A0_invdyn_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRR1G2P2A0_invdyn_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRR1G2P2A0_invdyn_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:24:47
% EndTime: 2020-03-09 21:24:49
% DurationCPUTime: 2.26s
% Computational Cost: add. (17778->341), mult. (13119->513), div. (2034->7), fcn. (8124->74), ass. (0->220)
t172 = pkin(7) + qJ(3,3);
t144 = sin(t172);
t188 = sin(qJ(3,3));
t93 = 0.1e1 / (pkin(1) * t144 + t188 * pkin(2));
t246 = t93 / 0.2e1;
t173 = pkin(7) + qJ(3,2);
t145 = sin(t173);
t189 = sin(qJ(3,2));
t94 = 0.1e1 / (pkin(1) * t145 + t189 * pkin(2));
t245 = t94 / 0.2e1;
t174 = pkin(7) + qJ(3,1);
t146 = sin(t174);
t190 = sin(qJ(3,1));
t95 = 0.1e1 / (pkin(1) * t146 + t190 * pkin(2));
t244 = t95 / 0.2e1;
t203 = m(2) + m(3);
t261 = m(1) * rSges(1,2);
t260 = m(2) * rSges(2,2);
t259 = m(3) * rSges(3,2);
t258 = pkin(1) * sin(pkin(7));
t181 = cos(pkin(7));
t257 = pkin(1) * t181;
t256 = pkin(2) * t181;
t202 = xDP(1);
t207 = 0.1e1 / pkin(3);
t231 = t207 / 0.2e1;
t223 = t202 * t231;
t175 = qJ(1,3) + pkin(7);
t182 = legFrame(3,2);
t127 = t182 + t175;
t128 = -t182 + t175;
t156 = qJ(1,3) + t182;
t157 = qJ(1,3) - t182;
t118 = qJ(3,3) + t127;
t119 = qJ(3,3) + t128;
t73 = sin(t118) + sin(t119);
t49 = t73 * pkin(3) + (sin(t127) + sin(t128)) * pkin(2) + (sin(t156) + sin(t157)) * pkin(1);
t219 = t49 * t93 * t223;
t201 = xDP(2);
t224 = t201 * t231;
t76 = -cos(t119) + cos(t118);
t52 = -t76 * pkin(3) + (cos(t128) - cos(t127)) * pkin(2) + (cos(t157) - cos(t156)) * pkin(1);
t46 = t52 * t93 * t224;
t200 = xDP(3);
t230 = t200 * t207;
t153 = qJ(1,3) + t172;
t136 = cos(t153);
t150 = cos(t175);
t192 = cos(qJ(1,3));
t70 = -t192 * pkin(1) - pkin(2) * t150 - pkin(3) * t136;
t252 = t70 * t93;
t64 = t230 * t252;
t25 = t46 + t64 - t219;
t232 = t202 / 0.2e1;
t233 = t201 / 0.2e1;
t249 = t136 * t93;
t40 = t200 * t249 + (t232 * t73 + t233 * t76) * t93;
t16 = t25 + t40;
t255 = t16 * t25;
t176 = qJ(1,2) + pkin(7);
t183 = legFrame(2,2);
t129 = t183 + t176;
t130 = -t183 + t176;
t158 = qJ(1,2) + t183;
t159 = qJ(1,2) - t183;
t120 = qJ(3,2) + t129;
t121 = qJ(3,2) + t130;
t74 = sin(t120) + sin(t121);
t50 = t74 * pkin(3) + (sin(t129) + sin(t130)) * pkin(2) + (sin(t158) + sin(t159)) * pkin(1);
t218 = t50 * t94 * t223;
t77 = -cos(t121) + cos(t120);
t53 = -t77 * pkin(3) + (cos(t130) - cos(t129)) * pkin(2) + (cos(t159) - cos(t158)) * pkin(1);
t47 = t53 * t94 * t224;
t154 = qJ(1,2) + t173;
t137 = cos(t154);
t151 = cos(t176);
t194 = cos(qJ(1,2));
t71 = -t194 * pkin(1) - pkin(2) * t151 - pkin(3) * t137;
t251 = t71 * t94;
t65 = t230 * t251;
t26 = t47 + t65 - t218;
t248 = t137 * t94;
t41 = t200 * t248 + (t232 * t74 + t233 * t77) * t94;
t17 = t26 + t41;
t254 = t17 * t26;
t177 = qJ(1,1) + pkin(7);
t184 = legFrame(1,2);
t131 = t184 + t177;
t132 = -t184 + t177;
t160 = qJ(1,1) + t184;
t161 = qJ(1,1) - t184;
t122 = qJ(3,1) + t131;
t123 = qJ(3,1) + t132;
t75 = sin(t122) + sin(t123);
t51 = t75 * pkin(3) + (sin(t131) + sin(t132)) * pkin(2) + (sin(t160) + sin(t161)) * pkin(1);
t217 = t51 * t95 * t223;
t78 = -cos(t123) + cos(t122);
t54 = -t78 * pkin(3) + (cos(t132) - cos(t131)) * pkin(2) + (cos(t161) - cos(t160)) * pkin(1);
t48 = t54 * t95 * t224;
t155 = qJ(1,1) + t174;
t138 = cos(t155);
t152 = cos(t177);
t196 = cos(qJ(1,1));
t72 = -t196 * pkin(1) - pkin(2) * t152 - pkin(3) * t138;
t250 = t72 * t95;
t66 = t230 * t250;
t27 = t48 + t66 - t217;
t247 = t138 * t95;
t42 = t200 * t247 + (t232 * t75 + t233 * t78) * t95;
t18 = t27 + t42;
t253 = t18 * t27;
t243 = t207 * t49;
t242 = t207 * t50;
t241 = t207 * t51;
t240 = t207 * t52;
t239 = t207 * t53;
t238 = t207 * t54;
t227 = pkin(1) * t259;
t113 = t144 * t227;
t226 = pkin(2) * t259;
t140 = t188 * t226;
t168 = rSges(3,1) ^ 2 + rSges(3,2) ^ 2;
t147 = cos(t172);
t191 = cos(qJ(3,3));
t222 = pkin(1) * t147 + pkin(2) * t191;
t213 = t222 * rSges(3,1);
t55 = Icges(3,3) - t113 - t140 + (t168 + t213) * m(3);
t237 = t207 * t55;
t114 = t145 * t227;
t141 = t189 * t226;
t148 = cos(t173);
t193 = cos(qJ(3,2));
t221 = pkin(1) * t148 + pkin(2) * t193;
t212 = t221 * rSges(3,1);
t56 = Icges(3,3) - t114 - t141 + (t168 + t212) * m(3);
t236 = t207 * t56;
t115 = t146 * t227;
t142 = t190 * t226;
t149 = cos(t174);
t195 = cos(qJ(3,1));
t220 = pkin(1) * t149 + pkin(2) * t195;
t211 = t220 * rSges(3,1);
t57 = Icges(3,3) - t115 - t142 + (t168 + t211) * m(3);
t235 = t207 * t57;
t117 = t168 * m(3) + Icges(3,3);
t234 = t117 * t207;
t229 = 0.2e1 * pkin(2) * pkin(3);
t209 = pkin(1) ^ 2;
t171 = pkin(2) ^ 2 + t209;
t228 = 0.2e1 * pkin(1);
t162 = sin(t182);
t163 = sin(t183);
t164 = sin(t184);
t165 = cos(t182);
t166 = cos(t183);
t167 = cos(t184);
t225 = (t162 * t165 + t163 * t166 + t164 * t167) * t203;
t139 = pkin(2) + t257;
t197 = rSges(3,2) * g(3);
t198 = rSges(3,1) * g(3);
t90 = t165 * g(1) - t162 * g(2);
t216 = sin(t153) * (rSges(3,2) * t90 + t198) - t136 * (rSges(3,1) * t90 - t197);
t91 = t166 * g(1) - t163 * g(2);
t215 = sin(t154) * (rSges(3,2) * t91 + t198) - t137 * (rSges(3,1) * t91 - t197);
t92 = t167 * g(1) - t164 * g(2);
t214 = sin(t155) * (rSges(3,2) * t92 + t198) - t138 * (rSges(3,1) * t92 - t197);
t143 = m(2) * rSges(2,1) + m(3) * pkin(2);
t210 = Icges(1,3) + Icges(2,3) + Icges(3,3) + (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) + (rSges(2,1) ^ 2 + t209 + (-0.2e1 * t258 + rSges(2,2)) * rSges(2,2)) * m(2) + 0.2e1 * t143 * t257;
t206 = pkin(3) ^ 2;
t187 = xDDP(1);
t186 = xDDP(2);
t185 = xDDP(3);
t179 = g(3) * t261;
t178 = g(3) * t260;
t125 = t143 * g(3);
t124 = t171 + t168;
t116 = m(1) * rSges(1,1) + t203 * pkin(1);
t100 = t116 * g(3);
t89 = t164 * g(1) + t167 * g(2);
t88 = t163 * g(1) + t166 * g(2);
t87 = t162 * g(1) + t165 * g(2);
t86 = t139 * rSges(3,1) - rSges(3,2) * t258;
t85 = rSges(3,1) * t258 + t139 * rSges(3,2);
t45 = -0.2e1 * t115 - 0.2e1 * t142 + (t124 + 0.2e1 * t211) * m(3) + t210;
t44 = -0.2e1 * t114 - 0.2e1 * t141 + (t124 + 0.2e1 * t212) * m(3) + t210;
t43 = -0.2e1 * t113 - 0.2e1 * t140 + (t124 + 0.2e1 * t213) * m(3) + t210;
t39 = (t138 * t57 + t72 * t234) * t95;
t38 = (t137 * t56 + t71 * t234) * t94;
t37 = (t136 * t55 + t70 * t234) * t93;
t36 = (t54 * t234 + t57 * t78) * t244;
t35 = (t53 * t234 + t56 * t77) * t245;
t34 = (t52 * t234 + t55 * t76) * t246;
t33 = (-t51 * t234 + t57 * t75) * t244;
t32 = (-t50 * t234 + t56 * t74) * t245;
t31 = (-t49 * t234 + t55 * t73) * t246;
t30 = (t138 * t45 + t72 * t235) * t95;
t29 = (t137 * t44 + t71 * t236) * t94;
t28 = (t136 * t43 + t70 * t237) * t93;
t24 = (t54 * t235 + t45 * t78) * t244;
t23 = (t53 * t236 + t44 * t77) * t245;
t22 = (t52 * t237 + t43 * t76) * t246;
t21 = (-t51 * t235 + t45 * t75) * t244;
t20 = (-t50 * t236 + t44 * t74) * t245;
t19 = (-t49 * t237 + t43 * t73) * t246;
t15 = -t217 / 0.2e1 + t48 / 0.2e1 + t66 / 0.2e1 + t42;
t14 = -t218 / 0.2e1 + t47 / 0.2e1 + t65 / 0.2e1 + t41;
t13 = -t219 / 0.2e1 + t46 / 0.2e1 + t64 / 0.2e1 + t40;
t12 = (-pkin(3) * t253 + (-pkin(3) * t18 - t220 * t42) * t42) * t95;
t11 = (-pkin(3) * t255 + (-pkin(3) * t16 - t222 * t40) * t40) * t93;
t10 = (-pkin(3) * t254 + (-pkin(3) * t17 - t221 * t41) * t41) * t94;
t9 = (t15 * t195 * t229 + t42 * t171 + t18 * t206 + (pkin(3) * t149 * t15 + t42 * t256) * t228) * t207 * t95 * t42 + (t139 * t195 - t190 * t258 + pkin(3)) / (t139 * t190 + t195 * t258) * t253;
t8 = (t14 * t193 * t229 + t17 * t206 + t41 * t171 + (pkin(3) * t14 * t148 + t41 * t256) * t228) * t207 * t94 * t41 + (t139 * t193 - t189 * t258 + pkin(3)) / (t139 * t189 + t193 * t258) * t254;
t7 = (t13 * t191 * t229 + t16 * t206 + t40 * t171 + (pkin(3) * t13 * t147 + t40 * t256) * t228) * t207 * t93 * t40 + (t139 * t191 - t188 * t258 + pkin(3)) / (t139 * t188 + t191 * t258) * t255;
t6 = -t117 * t9 - t57 * t12 + ((pkin(2) * (t190 * rSges(3,1) + rSges(3,2) * t195) + (rSges(3,1) * t146 + rSges(3,2) * t149) * pkin(1)) * t42 ^ 2 + t214) * m(3);
t5 = -t56 * t10 - t117 * t8 + ((pkin(2) * (t189 * rSges(3,1) + rSges(3,2) * t193) + (rSges(3,1) * t145 + rSges(3,2) * t148) * pkin(1)) * t41 ^ 2 + t215) * m(3);
t4 = -t55 * t11 - t117 * t7 + ((pkin(2) * (t188 * rSges(3,1) + rSges(3,2) * t191) + (rSges(3,1) * t144 + rSges(3,2) * t147) * pkin(1)) * t40 ^ 2 + t216) * m(3);
t3 = -t45 * t12 - t57 * t9 + (-t92 * t143 + t178) * t152 + (t92 * t260 + t125) * sin(t177) + (-t92 * t116 + t179) * t196 + (t92 * t261 + t100) * sin(qJ(1,1)) + (-0.2e1 * t27 * (t190 * t86 + t85 * t195) * t15 + t214) * m(3);
t2 = -t44 * t10 - t56 * t8 + (-t91 * t143 + t178) * t151 + (t91 * t260 + t125) * sin(t176) + (-t91 * t116 + t179) * t194 + (t91 * t261 + t100) * sin(qJ(1,2)) + (-0.2e1 * t26 * (t189 * t86 + t85 * t193) * t14 + t215) * m(3);
t1 = -t43 * t11 - t55 * t7 + (-t90 * t143 + t178) * t150 + (t90 * t260 + t125) * sin(t175) + (-t90 * t116 + t179) * t192 + (t90 * t261 + t100) * sin(qJ(1,3)) + (-0.2e1 * t25 * (t188 * t86 + t85 * t191) * t13 + t216) * m(3);
t58 = [(-g(1) + t187) * m(4) + t225 * t186 + (t19 * t249 + t20 * t248 + t21 * t247 + (t33 * t250 + t32 * t251 + t31 * t252) * t207) * t185 + (-t162 * t87 - t163 * t88 - t164 * t89 + (t162 ^ 2 + t163 ^ 2 + t164 ^ 2) * t187) * t203 + ((t21 * t75 - t33 * t241) * t187 + (t21 * t78 + t33 * t238) * t186 + t75 * t3 - t6 * t241) * t244 + ((t20 * t74 - t32 * t242) * t187 + (t20 * t77 + t32 * t239) * t186 + t74 * t2 - t5 * t242) * t245 + ((t19 * t73 - t31 * t243) * t187 + (t19 * t76 + t31 * t240) * t186 + t73 * t1 - t4 * t243) * t246; (-g(2) + t186) * m(4) + t225 * t187 + (t22 * t249 + t23 * t248 + t24 * t247 + (t36 * t250 + t35 * t251 + t34 * t252) * t207) * t185 + (-t165 * t87 - t166 * t88 - t167 * t89 + (t165 ^ 2 + t166 ^ 2 + t167 ^ 2) * t186) * t203 + ((t24 * t75 - t36 * t241) * t187 + (t36 * t238 + t24 * t78) * t186 + t78 * t3 + t6 * t238) * t244 + ((t23 * t74 - t35 * t242) * t187 + (t23 * t77 + t35 * t239) * t186 + t77 * t2 + t5 * t239) * t245 + ((t22 * t73 - t34 * t243) * t187 + (t22 * t76 + t34 * t240) * t186 + t76 * t1 + t4 * t240) * t246; t1 * t249 + t2 * t248 + t3 * t247 - m(4) * g(3) + (t6 * t250 + t5 * t251 + t4 * t252) * t207 + (t28 * t249 + t29 * t248 + t30 * t247 + m(4) + (t39 * t250 + t38 * t251 + t37 * t252) * t207) * t185 + ((-t39 * t241 + t30 * t75) * t187 + (t39 * t238 + t30 * t78) * t186) * t244 + ((-t38 * t242 + t29 * t74) * t187 + (t38 * t239 + t29 * t77) * t186) * t245 + ((-t37 * t243 + t28 * t73) * t187 + (t37 * t240 + t28 * t76) * t186) * t246;];
tauX  = t58;
