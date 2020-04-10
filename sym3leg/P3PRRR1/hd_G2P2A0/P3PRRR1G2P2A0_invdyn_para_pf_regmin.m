% Calculate minimal parameter regressor of inverse dynamics forces for
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tauX_reg [3x8]
%   minimal parameter regressor of inverse dynamics force vector
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:18
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX_reg = P3PRRR1G2P2A0_invdyn_para_pf_regmin(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR1G2P2A0_invdyn_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRR1G2P2A0_invdyn_para_pf_regmin: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3PRRR1G2P2A0_invdyn_para_pf_regmin: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR1G2P2A0_invdyn_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRR1G2P2A0_invdyn_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3PRRR1G2P2A0_invdyn_para_pf_regmin: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR1G2P2A0_invdyn_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR1G2P2A0_invdyn_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_reg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:18:25
% EndTime: 2020-03-09 21:18:27
% DurationCPUTime: 2.08s
% Computational Cost: add. (23254->238), mult. (13422->448), div. (3195->10), fcn. (13746->36), ass. (0->219)
t171 = pkin(7) + qJ(2,1);
t162 = qJ(3,1) + t171;
t153 = cos(t162);
t159 = cos(t171);
t150 = sin(t162);
t156 = sin(t171);
t81 = -t159 * t150 + t153 * t156;
t282 = 0.1e1 / t81;
t256 = t282 * (pkin(2) * t159 + pkin(3) * t153);
t170 = pkin(7) + qJ(2,2);
t161 = qJ(3,2) + t170;
t152 = cos(t161);
t158 = cos(t170);
t149 = sin(t161);
t155 = sin(t170);
t80 = -t158 * t149 + t152 * t155;
t283 = 0.1e1 / t80;
t257 = t283 * (pkin(2) * t158 + pkin(3) * t152);
t169 = pkin(7) + qJ(2,3);
t160 = qJ(3,3) + t169;
t151 = cos(t160);
t157 = cos(t169);
t148 = sin(t160);
t154 = sin(t169);
t79 = -t157 * t148 + t151 * t154;
t284 = 0.1e1 / t79;
t258 = t284 * (pkin(2) * t157 + pkin(3) * t151);
t172 = legFrame(3,2);
t163 = sin(t172);
t250 = t163 * t284;
t173 = legFrame(2,2);
t164 = sin(t173);
t249 = t164 * t283;
t174 = legFrame(1,2);
t165 = sin(t174);
t248 = t165 * t282;
t166 = cos(t172);
t287 = t166 * t284;
t167 = cos(t173);
t286 = t167 * t283;
t168 = cos(t174);
t285 = t168 * t282;
t188 = 0.1e1 / pkin(3);
t281 = 0.2e1 * t188;
t280 = -g(1) / 0.2e1;
t279 = g(1) / 0.2e1;
t278 = -g(2) / 0.2e1;
t277 = g(2) / 0.2e1;
t191 = 0.1e1 / pkin(2) ^ 2;
t184 = xDP(3);
t185 = xDP(2);
t186 = xDP(1);
t85 = -t163 * t185 + t166 * t186;
t200 = -t85 * t148 - t151 * t184;
t219 = t148 * t287;
t177 = xDDP(1);
t190 = 0.1e1 / pkin(2);
t236 = t177 * t190;
t43 = t219 * t236;
t176 = xDDP(2);
t237 = t176 * t190;
t212 = t237 * t250;
t46 = t148 * t212;
t175 = xDDP(3);
t238 = t175 * t190;
t253 = t151 * t284;
t55 = t238 * t253;
t232 = -t43 + t46 - t55;
t197 = t148 * t154 + t151 * t157;
t264 = t200 * t284;
t34 = -pkin(2) * (t85 * t154 + t157 * t184) + t200 * pkin(3);
t241 = t34 * t188;
t244 = t190 * t284;
t31 = (-t200 + t241) * t244;
t25 = -t31 * pkin(3) + t197 * t264;
t261 = t284 ^ 2;
t267 = t31 * t34;
t19 = -(t200 * t25 + t267) * t191 * t261 + t232;
t276 = pkin(2) * t19;
t86 = -t164 * t185 + t167 * t186;
t199 = -t86 * t149 - t152 * t184;
t217 = t149 * t286;
t44 = t217 * t236;
t211 = t237 * t249;
t47 = t149 * t211;
t252 = t152 * t283;
t56 = t238 * t252;
t231 = -t44 + t47 - t56;
t196 = t149 * t155 + t152 * t158;
t263 = t199 * t283;
t35 = -pkin(2) * (t86 * t155 + t158 * t184) + t199 * pkin(3);
t240 = t35 * t188;
t243 = t190 * t283;
t32 = (-t199 + t240) * t243;
t26 = -t32 * pkin(3) + t196 * t263;
t260 = t283 ^ 2;
t266 = t32 * t35;
t20 = -(t199 * t26 + t266) * t191 * t260 + t231;
t275 = pkin(2) * t20;
t87 = -t165 * t185 + t168 * t186;
t198 = -t87 * t150 - t153 * t184;
t215 = t150 * t285;
t45 = t215 * t236;
t210 = t237 * t248;
t48 = t150 * t210;
t251 = t153 * t282;
t57 = t238 * t251;
t230 = -t45 + t48 - t57;
t259 = t282 ^ 2;
t36 = -pkin(2) * (t87 * t156 + t159 * t184) + t198 * pkin(3);
t239 = t36 * t188;
t242 = t190 * t282;
t33 = (-t198 + t239) * t242;
t265 = t33 * t36;
t195 = t150 * t156 + t153 * t159;
t262 = t198 * t282;
t27 = -t33 * pkin(3) + t195 * t262;
t21 = -(t198 * t27 + t265) * t191 * t259 + t230;
t274 = pkin(2) * t21;
t142 = t172 + t160;
t273 = sin(t142) / 0.2e1;
t144 = t173 + t161;
t272 = sin(t144) / 0.2e1;
t146 = t174 + t162;
t271 = sin(t146) / 0.2e1;
t143 = -t172 + t160;
t270 = -cos(t143) / 0.2e1;
t145 = -t173 + t161;
t269 = -cos(t145) / 0.2e1;
t147 = -t174 + t162;
t268 = -cos(t147) / 0.2e1;
t255 = t176 - g(2);
t254 = t177 - g(1);
t187 = pkin(3) ^ 2;
t234 = 0.2e1 * pkin(3);
t28 = (-t200 + t241 / 0.2e1) * t244;
t247 = t188 * (t31 * t187 + (t197 * t28 * t234 - t264) * pkin(2));
t29 = (-t199 + t240 / 0.2e1) * t243;
t246 = t188 * (t32 * t187 + (t196 * t29 * t234 - t263) * pkin(2));
t30 = (-t198 + t239 / 0.2e1) * t242;
t245 = t188 * (t33 * t187 + (t195 * t30 * t234 - t262) * pkin(2));
t235 = t188 * t190;
t233 = 0.2e1 * t235;
t88 = pkin(2) * t154 + pkin(3) * t148;
t229 = t88 * t250;
t89 = pkin(2) * t155 + pkin(3) * t149;
t228 = t89 * t249;
t90 = pkin(2) * t156 + pkin(3) * t150;
t227 = t90 * t248;
t226 = t88 * t287;
t225 = t89 * t286;
t224 = t90 * t285;
t74 = 0.1e1 / t79 ^ 2;
t223 = t190 * t200 ^ 2 * t74;
t76 = 0.1e1 / t80 ^ 2;
t222 = t190 * t199 ^ 2 * t76;
t78 = 0.1e1 / t81 ^ 2;
t221 = t190 * t198 ^ 2 * t78;
t220 = t148 * t250;
t218 = t149 * t249;
t216 = t150 * t248;
t214 = t175 * t235;
t213 = t177 * t235;
t209 = t31 * (t197 * pkin(2) + pkin(3)) * t74 * t241;
t208 = t32 * (t196 * pkin(2) + pkin(3)) * t76 * t240;
t207 = t33 * (t195 * pkin(2) + pkin(3)) * t78 * t239;
t125 = sin(t143);
t130 = cos(t142);
t206 = g(1) * t270 + g(2) * t273 + g(3) * t148 + t125 * t278 + t130 * t280;
t127 = sin(t145);
t132 = cos(t144);
t205 = g(1) * t269 + g(2) * t272 + g(3) * t149 + t127 * t278 + t132 * t280;
t129 = sin(t147);
t134 = cos(t146);
t204 = g(1) * t268 + g(2) * t271 + g(3) * t150 + t129 * t278 + t134 * t280;
t203 = g(1) * t273 + g(2) * t270 + g(3) * t151 + t125 * t279 + t130 * t277;
t202 = g(1) * t272 + g(2) * t269 + g(3) * t152 + t127 * t279 + t132 * t277;
t201 = g(1) * t271 + g(2) * t268 + g(3) * t153 + t129 * t279 + t134 * t277;
t194 = -t188 * t88 * t212 + t213 * t226 + t214 * t258;
t193 = -t188 * t89 * t211 + t213 * t225 + t214 * t257;
t192 = -t188 * t90 * t210 + t213 * t224 + t214 * t256;
t189 = 0.1e1 / pkin(3) ^ 2;
t183 = cos(qJ(3,1));
t182 = cos(qJ(3,2));
t181 = cos(qJ(3,3));
t180 = sin(qJ(3,1));
t179 = sin(qJ(3,2));
t178 = sin(qJ(3,3));
t96 = -t168 * g(1) + t165 * g(2);
t95 = -t167 * g(1) + t164 * g(2);
t94 = -t166 * g(1) + t163 * g(2);
t69 = g(3) * t159 - t96 * t156;
t68 = g(3) * t158 - t95 * t155;
t67 = g(3) * t157 - t94 * t154;
t66 = g(3) * t156 + t96 * t159;
t65 = g(3) * t155 + t95 * t158;
t64 = g(3) * t154 + t94 * t157;
t60 = t254 * t165 + t255 * t168;
t59 = t254 * t164 + t255 * t167;
t58 = t254 * t163 + t255 * t166;
t18 = t180 * t221 + t183 * t274 + t204;
t17 = t179 * t222 + t182 * t275 + t205;
t16 = t178 * t223 + t181 * t276 + t206;
t15 = -t180 * t274 + t183 * t221 + t201;
t14 = -t179 * t275 + t182 * t222 + t202;
t13 = -t178 * t276 + t181 * t223 + t203;
t12 = (t207 - (t265 - (-t27 - t245) * t198) * t259) * t191 + t192 + t230;
t11 = -0.2e1 * t45 + 0.2e1 * t48 - 0.2e1 * t57 + (t207 - (0.2e1 * t265 - (-0.2e1 * t27 - t245) * t198) * t259) * t191 + t192;
t10 = (t208 - (t266 - (-t26 - t246) * t199) * t260) * t191 + t193 + t231;
t9 = -0.2e1 * t44 + 0.2e1 * t47 - 0.2e1 * t56 + (t208 - (0.2e1 * t266 - (-0.2e1 * t26 - t246) * t199) * t260) * t191 + t193;
t8 = (t209 - (t267 - (-t25 - t247) * t200) * t261) * t191 + t194 + t232;
t7 = -0.2e1 * t43 + 0.2e1 * t46 - 0.2e1 * t55 + (t209 - (0.2e1 * t267 - (-0.2e1 * t25 - t247) * t200) * t261) * t191 + t194;
t6 = -pkin(2) * (t180 * t11 + (t189 * t36 - t198 * t281) * t78 * t36 * t183 * t191) + t201;
t5 = -pkin(2) * (t179 * t9 + (t189 * t35 - t199 * t281) * t76 * t35 * t182 * t191) + t202;
t4 = -pkin(2) * (t178 * t7 + (t189 * t34 - t200 * t281) * t74 * t34 * t181 * t191) + t203;
t3 = pkin(2) * (-t180 * t233 * t282 * t30 * t36 + t11 * t183) + t204;
t2 = pkin(2) * (-t179 * t233 * t283 * t29 * t35 + t9 * t182) + t205;
t1 = pkin(2) * (-t178 * t233 * t28 * t284 * t34 + t7 * t181) + t206;
t22 = [t163 * t58 + t164 * t59 + t165 * t60, (-t19 * t219 - t20 * t217 - t21 * t215) * t190, (-t66 * t215 - t65 * t217 - t64 * t219) * t190, (-t69 * t215 - t68 * t217 - t67 * t219) * t190, (-t10 * t217 - t12 * t215 - t8 * t219 + (t10 * t225 + t12 * t224 + t8 * t226) * t188) * t190, (-t1 * t219 - t2 * t217 - t3 * t215 + (t16 * t226 + t17 * t225 + t18 * t224) * t188) * t190, (-t4 * t219 - t5 * t217 - t6 * t215 + (t13 * t226 + t14 * t225 + t15 * t224) * t188) * t190, t254; t166 * t58 + t167 * t59 + t168 * t60, (t19 * t220 + t20 * t218 + t21 * t216) * t190, (t66 * t216 + t65 * t218 + t64 * t220) * t190, (t69 * t216 + t68 * t218 + t67 * t220) * t190, (t10 * t218 + t12 * t216 + t8 * t220 + (-t10 * t228 - t12 * t227 - t8 * t229) * t188) * t190, (t1 * t220 + t2 * t218 + t3 * t216 + (-t16 * t229 - t17 * t228 - t18 * t227) * t188) * t190, (t4 * t220 + t5 * t218 + t6 * t216 + (-t13 * t229 - t14 * t228 - t15 * t227) * t188) * t190, t255; 0, (-t19 * t253 - t20 * t252 - t21 * t251) * t190, (-t66 * t251 - t65 * t252 - t64 * t253) * t190, (-t69 * t251 - t68 * t252 - t67 * t253) * t190, (-t10 * t252 - t12 * t251 - t8 * t253 + (t10 * t257 + t12 * t256 + t8 * t258) * t188) * t190, (-t1 * t253 - t2 * t252 - t3 * t251 + (t16 * t258 + t17 * t257 + t18 * t256) * t188) * t190, (-t4 * t253 - t5 * t252 - t6 * t251 + (t13 * t258 + t14 * t257 + t15 * t256) * t188) * t190, t175 - g(3);];
tauX_reg  = t22;
