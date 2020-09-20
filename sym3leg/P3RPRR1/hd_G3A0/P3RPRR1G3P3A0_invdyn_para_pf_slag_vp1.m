% Calculate vector of inverse dynamics forces for parallel robot
% P3RPRR1G3P3A0
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
% Datum: 2020-03-09 21:27
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3RPRR1G3P3A0_invdyn_para_pf_slag_vp1(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRR1G3P3A0_invdyn_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRR1G3P3A0_invdyn_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RPRR1G3P3A0_invdyn_para_pf_slag_vp1: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRR1G3P3A0_invdyn_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRR1G3P3A0_invdyn_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRR1G3P3A0_invdyn_para_pf_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRR1G3P3A0_invdyn_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RPRR1G3P3A0_invdyn_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RPRR1G3P3A0_invdyn_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRR1G3P3A0_invdyn_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRR1G3P3A0_invdyn_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:26:43
% EndTime: 2020-03-09 21:26:46
% DurationCPUTime: 2.36s
% Computational Cost: add. (17778->353), mult. (13119->526), div. (2034->7), fcn. (8124->74), ass. (0->230)
t169 = pkin(7) + qJ(3,3);
t141 = sin(t169);
t183 = sin(qJ(3,3));
t90 = 0.1e1 / (pkin(1) * t141 + t183 * pkin(2));
t245 = t90 / 0.2e1;
t170 = pkin(7) + qJ(3,2);
t142 = sin(t170);
t185 = sin(qJ(3,2));
t91 = 0.1e1 / (pkin(1) * t142 + t185 * pkin(2));
t244 = t91 / 0.2e1;
t171 = pkin(7) + qJ(3,1);
t143 = sin(t171);
t187 = sin(qJ(3,1));
t92 = 0.1e1 / (pkin(1) * t143 + t187 * pkin(2));
t243 = t92 / 0.2e1;
t196 = m(2) + m(3);
t268 = m(1) * rSges(1,2);
t267 = m(2) * rSges(2,2);
t266 = m(3) * rSges(3,2);
t265 = rSges(3,1) * g(3);
t264 = rSges(3,2) * g(3);
t263 = pkin(1) * sin(pkin(7));
t176 = cos(pkin(7));
t262 = pkin(1) * t176;
t261 = pkin(2) * t176;
t194 = xDP(2);
t200 = 0.1e1 / pkin(3);
t227 = t200 / 0.2e1;
t218 = t194 * t227;
t172 = qJ(1,3) + pkin(7);
t177 = legFrame(3,2);
t124 = t177 + t172;
t125 = -t177 + t172;
t153 = qJ(1,3) + t177;
t154 = qJ(1,3) - t177;
t115 = qJ(3,3) + t124;
t116 = qJ(3,3) + t125;
t70 = -sin(t115) + sin(t116);
t49 = t70 * pkin(3) + (-sin(t124) + sin(t125)) * pkin(2) + (-sin(t153) + sin(t154)) * pkin(1);
t216 = t49 * t218;
t150 = qJ(1,3) + t169;
t130 = sin(t150);
t193 = xDP(3);
t232 = t130 * t193;
t195 = xDP(1);
t228 = t195 / 0.2e1;
t229 = t194 / 0.2e1;
t73 = cos(t116) + cos(t115);
t251 = (t228 * t73 + t229 * t70) * t90;
t217 = t195 * t227;
t52 = -t73 * pkin(3) + (-cos(t125) - cos(t124)) * pkin(2) + (-cos(t154) - cos(t153)) * pkin(1);
t46 = t52 * t90 * t217;
t226 = t193 * t200;
t144 = sin(t172);
t184 = sin(qJ(1,3));
t67 = t184 * pkin(1) + pkin(2) * t144 + pkin(3) * t130;
t257 = t67 * t90;
t64 = t226 * t257;
t254 = t46 + t64;
t16 = (-t216 - t232) * t90 + t251 + t254;
t25 = -t90 * t216 + t254;
t260 = t16 * t25;
t173 = qJ(1,2) + pkin(7);
t178 = legFrame(2,2);
t126 = t178 + t173;
t127 = -t178 + t173;
t155 = qJ(1,2) + t178;
t156 = qJ(1,2) - t178;
t117 = qJ(3,2) + t126;
t118 = qJ(3,2) + t127;
t71 = -sin(t117) + sin(t118);
t50 = t71 * pkin(3) + (-sin(t126) + sin(t127)) * pkin(2) + (-sin(t155) + sin(t156)) * pkin(1);
t215 = t50 * t218;
t151 = qJ(1,2) + t170;
t131 = sin(t151);
t231 = t131 * t193;
t74 = cos(t118) + cos(t117);
t250 = (t228 * t74 + t229 * t71) * t91;
t53 = -t74 * pkin(3) + (-cos(t127) - cos(t126)) * pkin(2) + (-cos(t156) - cos(t155)) * pkin(1);
t47 = t53 * t91 * t217;
t145 = sin(t173);
t186 = sin(qJ(1,2));
t68 = t186 * pkin(1) + pkin(2) * t145 + pkin(3) * t131;
t256 = t68 * t91;
t65 = t226 * t256;
t253 = t47 + t65;
t17 = (-t215 - t231) * t91 + t250 + t253;
t26 = -t91 * t215 + t253;
t259 = t17 * t26;
t174 = qJ(1,1) + pkin(7);
t179 = legFrame(1,2);
t128 = t179 + t174;
t129 = -t179 + t174;
t157 = qJ(1,1) + t179;
t158 = qJ(1,1) - t179;
t119 = qJ(3,1) + t128;
t120 = qJ(3,1) + t129;
t72 = -sin(t119) + sin(t120);
t51 = t72 * pkin(3) + (-sin(t128) + sin(t129)) * pkin(2) + (-sin(t157) + sin(t158)) * pkin(1);
t214 = t51 * t218;
t152 = qJ(1,1) + t171;
t132 = sin(t152);
t230 = t132 * t193;
t75 = cos(t120) + cos(t119);
t249 = (t228 * t75 + t229 * t72) * t92;
t54 = -t75 * pkin(3) + (-cos(t129) - cos(t128)) * pkin(2) + (-cos(t158) - cos(t157)) * pkin(1);
t48 = t54 * t92 * t217;
t146 = sin(t174);
t188 = sin(qJ(1,1));
t69 = t188 * pkin(1) + pkin(2) * t146 + pkin(3) * t132;
t255 = t69 * t92;
t66 = t226 * t255;
t252 = t48 + t66;
t18 = (-t214 - t230) * t92 + t249 + t252;
t27 = -t92 * t214 + t252;
t258 = t18 * t27;
t248 = t130 * t90;
t247 = t131 * t91;
t246 = t132 * t92;
t242 = t200 * t49;
t241 = t200 * t50;
t240 = t200 * t51;
t239 = t200 * t52;
t238 = t200 * t53;
t237 = t200 * t54;
t222 = pkin(1) * t266;
t110 = t141 * t222;
t221 = pkin(2) * t266;
t137 = t183 * t221;
t165 = rSges(3,1) ^ 2 + rSges(3,2) ^ 2;
t147 = cos(t169);
t189 = cos(qJ(3,3));
t212 = pkin(1) * t147 + pkin(2) * t189;
t206 = t212 * rSges(3,1);
t55 = Icges(3,3) - t110 - t137 + (t165 + t206) * m(3);
t236 = t200 * t55;
t111 = t142 * t222;
t138 = t185 * t221;
t148 = cos(t170);
t190 = cos(qJ(3,2));
t211 = pkin(1) * t148 + pkin(2) * t190;
t205 = t211 * rSges(3,1);
t56 = Icges(3,3) - t111 - t138 + (t165 + t205) * m(3);
t235 = t200 * t56;
t112 = t143 * t222;
t139 = t187 * t221;
t149 = cos(t171);
t191 = cos(qJ(3,1));
t210 = pkin(1) * t149 + pkin(2) * t191;
t204 = t210 * rSges(3,1);
t57 = Icges(3,3) - t112 - t139 + (t165 + t204) * m(3);
t234 = t200 * t57;
t114 = t165 * m(3) + Icges(3,3);
t233 = t114 * t200;
t225 = 0.2e1 * pkin(2) * pkin(3);
t202 = pkin(1) ^ 2;
t168 = pkin(2) ^ 2 + t202;
t224 = 0.2e1 * pkin(1);
t223 = g(3) * t268;
t220 = g(3) * t267;
t159 = sin(t177);
t160 = sin(t178);
t161 = sin(t179);
t162 = cos(t177);
t163 = cos(t178);
t164 = cos(t179);
t219 = (t159 * t162 + t160 * t163 + t161 * t164) * t196;
t136 = pkin(2) + t262;
t213 = -t218 / 0.2e1;
t87 = t162 * g(1) - t159 * g(2);
t209 = t130 * (rSges(3,1) * t87 - t264) + cos(t150) * (rSges(3,2) * t87 + t265);
t88 = t163 * g(1) - t160 * g(2);
t208 = t131 * (rSges(3,1) * t88 - t264) + cos(t151) * (rSges(3,2) * t88 + t265);
t89 = t164 * g(1) - t161 * g(2);
t207 = t132 * (rSges(3,1) * t89 - t264) + cos(t152) * (rSges(3,2) * t89 + t265);
t140 = m(2) * rSges(2,1) + m(3) * pkin(2);
t203 = Icges(1,3) + Icges(2,3) + Icges(3,3) + (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) + (rSges(2,1) ^ 2 + t202 + (-0.2e1 * t263 + rSges(2,2)) * rSges(2,2)) * m(2) + 0.2e1 * t140 * t262;
t199 = pkin(3) ^ 2;
t182 = xDDP(1);
t181 = xDDP(2);
t180 = xDDP(3);
t122 = g(3) * t140;
t121 = t168 + t165;
t113 = m(1) * rSges(1,1) + t196 * pkin(1);
t97 = g(3) * t113;
t86 = t161 * g(1) + t164 * g(2);
t85 = t160 * g(1) + t163 * g(2);
t84 = t159 * g(1) + t162 * g(2);
t83 = t136 * rSges(3,1) - rSges(3,2) * t263;
t82 = rSges(3,1) * t263 + t136 * rSges(3,2);
t45 = -0.2e1 * t112 - 0.2e1 * t139 + (t121 + 0.2e1 * t204) * m(3) + t203;
t44 = -0.2e1 * t111 - 0.2e1 * t138 + (t121 + 0.2e1 * t205) * m(3) + t203;
t43 = -0.2e1 * t110 - 0.2e1 * t137 + (t121 + 0.2e1 * t206) * m(3) + t203;
t42 = -t92 * t230 + t249;
t41 = -t91 * t231 + t250;
t40 = -t90 * t232 + t251;
t39 = (-t132 * t57 + t69 * t233) * t92;
t38 = (-t131 * t56 + t68 * t233) * t91;
t37 = (-t130 * t55 + t67 * t233) * t90;
t36 = (t54 * t233 + t57 * t75) * t243;
t35 = (t53 * t233 + t56 * t74) * t244;
t34 = (t52 * t233 + t55 * t73) * t245;
t33 = (-t51 * t233 + t57 * t72) * t243;
t32 = (-t50 * t233 + t56 * t71) * t244;
t31 = (-t49 * t233 + t55 * t70) * t245;
t30 = (-t132 * t45 + t69 * t234) * t92;
t29 = (-t131 * t44 + t68 * t235) * t91;
t28 = (-t130 * t43 + t67 * t236) * t90;
t24 = (t54 * t234 + t45 * t75) * t243;
t23 = (t53 * t235 + t44 * t74) * t244;
t22 = (t52 * t236 + t43 * t73) * t245;
t21 = (-t51 * t234 + t45 * t72) * t243;
t20 = (-t50 * t235 + t44 * t71) * t244;
t19 = (-t49 * t236 + t43 * t70) * t245;
t15 = t48 / 0.2e1 + t66 / 0.2e1 + (t51 * t213 - t230) * t92 + t249;
t14 = t47 / 0.2e1 + t65 / 0.2e1 + (t50 * t213 - t231) * t91 + t250;
t13 = t46 / 0.2e1 + t64 / 0.2e1 + (t49 * t213 - t232) * t90 + t251;
t12 = (-pkin(3) * t260 + (-pkin(3) * t16 - t212 * t40) * t40) * t90;
t11 = (-pkin(3) * t258 + (-pkin(3) * t18 - t210 * t42) * t42) * t92;
t10 = (-pkin(3) * t259 + (-pkin(3) * t17 - t211 * t41) * t41) * t91;
t9 = (t15 * t191 * t225 + t42 * t168 + t18 * t199 + (pkin(3) * t149 * t15 + t42 * t261) * t224) * t200 * t92 * t42 + (t136 * t191 - t187 * t263 + pkin(3)) / (t136 * t187 + t191 * t263) * t258;
t8 = (t14 * t190 * t225 + t41 * t168 + t17 * t199 + (pkin(3) * t14 * t148 + t41 * t261) * t224) * t200 * t91 * t41 + (t136 * t190 - t185 * t263 + pkin(3)) / (t136 * t185 + t190 * t263) * t259;
t7 = (t13 * t189 * t225 + t16 * t199 + t40 * t168 + (pkin(3) * t13 * t147 + t40 * t261) * t224) * t200 * t90 * t40 + (t136 * t189 - t183 * t263 + pkin(3)) / (t136 * t183 + t189 * t263) * t260;
t6 = -t57 * t11 - t114 * t9 + ((pkin(2) * (t187 * rSges(3,1) + rSges(3,2) * t191) + (rSges(3,1) * t143 + rSges(3,2) * t149) * pkin(1)) * t42 ^ 2 + t207) * m(3);
t5 = -t56 * t10 - t114 * t8 + ((pkin(2) * (t185 * rSges(3,1) + rSges(3,2) * t190) + (rSges(3,1) * t142 + rSges(3,2) * t148) * pkin(1)) * t41 ^ 2 + t208) * m(3);
t4 = -t114 * t7 - t55 * t12 + ((pkin(2) * (t183 * rSges(3,1) + rSges(3,2) * t189) + (rSges(3,1) * t141 + rSges(3,2) * t147) * pkin(1)) * t40 ^ 2 + t209) * m(3);
t3 = -t45 * t11 - t57 * t9 + (t89 * t267 + t122) * cos(t174) + (t140 * t89 - t220) * t146 + (t89 * t268 + t97) * cos(qJ(1,1)) + (t113 * t89 - t223) * t188 + (-0.2e1 * (t187 * t83 + t82 * t191) * t27 * t15 + t207) * m(3);
t2 = -t44 * t10 - t56 * t8 + (t88 * t267 + t122) * cos(t173) + (t140 * t88 - t220) * t145 + (t88 * t268 + t97) * cos(qJ(1,2)) + (t113 * t88 - t223) * t186 + (-0.2e1 * (t185 * t83 + t82 * t190) * t26 * t14 + t208) * m(3);
t1 = -t43 * t12 - t55 * t7 + (t87 * t267 + t122) * cos(t172) + (t140 * t87 - t220) * t144 + (t87 * t268 + t97) * cos(qJ(1,3)) + (t113 * t87 - t223) * t184 + (-0.2e1 * (t183 * t83 + t82 * t189) * t25 * t13 + t209) * m(3);
t58 = [(-g(1) + t182) * m(4) + t219 * t181 + (-t22 * t248 - t23 * t247 - t24 * t246 + (t36 * t255 + t35 * t256 + t34 * t257) * t200) * t180 + (-t159 * t84 - t160 * t85 - t161 * t86 + (t159 ^ 2 + t160 ^ 2 + t161 ^ 2) * t182) * t196 + ((t36 * t237 + t24 * t75) * t182 + (t24 * t72 - t36 * t240) * t181 + t75 * t3 + t6 * t237) * t243 + ((t23 * t74 + t35 * t238) * t182 + (t23 * t71 - t35 * t241) * t181 + t74 * t2 + t5 * t238) * t244 + ((t22 * t73 + t34 * t239) * t182 + (t22 * t70 - t34 * t242) * t181 + t73 * t1 + t4 * t239) * t245; (-g(2) + t181) * m(4) + t219 * t182 + (-t19 * t248 - t20 * t247 - t21 * t246 + (t33 * t255 + t32 * t256 + t31 * t257) * t200) * t180 + (-t162 * t84 - t163 * t85 - t164 * t86 + (t162 ^ 2 + t163 ^ 2 + t164 ^ 2) * t181) * t196 + ((t21 * t75 + t33 * t237) * t182 + (t21 * t72 - t33 * t240) * t181 + t72 * t3 - t6 * t240) * t243 + ((t20 * t74 + t32 * t238) * t182 + (t20 * t71 - t32 * t241) * t181 + t71 * t2 - t5 * t241) * t244 + ((t19 * t73 + t31 * t239) * t182 + (t19 * t70 - t31 * t242) * t181 + t70 * t1 - t4 * t242) * t245; -t1 * t248 - t2 * t247 - t3 * t246 - m(4) * g(3) + (t6 * t255 + t5 * t256 + t4 * t257) * t200 + (-t28 * t248 - t29 * t247 - t30 * t246 + m(4) + (t39 * t255 + t38 * t256 + t37 * t257) * t200) * t180 + ((t39 * t237 + t30 * t75) * t182 + (-t39 * t240 + t30 * t72) * t181) * t243 + ((t38 * t238 + t29 * t74) * t182 + (-t38 * t241 + t29 * t71) * t181) * t244 + ((t37 * t239 + t28 * t73) * t182 + (-t37 * t242 + t28 * t70) * t181) * t245;];
tauX  = t58;
