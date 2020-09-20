% Calculate minimal parameter regressor of inverse dynamics forces for
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
% Datum: 2020-03-09 21:07
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX_reg = P3PRRR1G3P3A0_invdyn_para_pf_regmin(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR1G3P3A0_invdyn_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRR1G3P3A0_invdyn_para_pf_regmin: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3PRRR1G3P3A0_invdyn_para_pf_regmin: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR1G3P3A0_invdyn_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRR1G3P3A0_invdyn_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3PRRR1G3P3A0_invdyn_para_pf_regmin: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR1G3P3A0_invdyn_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR1G3P3A0_invdyn_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_reg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:06:53
% EndTime: 2020-03-09 21:06:55
% DurationCPUTime: 2.31s
% Computational Cost: add. (22876->238), mult. (13368->438), div. (3141->10), fcn. (13530->36), ass. (0->212)
t163 = pkin(7) + qJ(2,3);
t154 = qJ(3,3) + t163;
t142 = sin(t154);
t145 = cos(t154);
t178 = xDP(3);
t166 = legFrame(3,2);
t157 = sin(t166);
t160 = cos(t166);
t179 = xDP(2);
t180 = xDP(1);
t200 = t157 * t179 - t160 * t180;
t55 = t142 * t178 + t200 * t145;
t164 = pkin(7) + qJ(2,2);
t155 = qJ(3,2) + t164;
t143 = sin(t155);
t146 = cos(t155);
t167 = legFrame(2,2);
t158 = sin(t167);
t161 = cos(t167);
t199 = t158 * t179 - t161 * t180;
t56 = t143 * t178 + t199 * t146;
t165 = pkin(7) + qJ(2,1);
t156 = qJ(3,1) + t165;
t144 = sin(t156);
t147 = cos(t156);
t168 = legFrame(1,2);
t159 = sin(t168);
t162 = cos(t168);
t198 = t159 * t179 - t162 * t180;
t57 = t144 * t178 + t198 * t147;
t184 = 1 / pkin(2);
t286 = -2 * t184;
t148 = sin(t163);
t151 = cos(t163);
t76 = t142 * t151 - t148 * t145;
t282 = 0.1e1 / t76;
t285 = t160 * t282;
t149 = sin(t164);
t152 = cos(t164);
t77 = t143 * t152 - t149 * t146;
t281 = 0.1e1 / t77;
t284 = t161 * t281;
t150 = sin(t165);
t153 = cos(t165);
t78 = t144 * t153 - t150 * t147;
t280 = 0.1e1 / t78;
t283 = t162 * t280;
t71 = 0.1e1 / t76 ^ 2;
t73 = 0.1e1 / t77 ^ 2;
t75 = 0.1e1 / t78 ^ 2;
t279 = g(1) / 0.2e1;
t278 = -g(2) / 0.2e1;
t185 = 0.1e1 / pkin(2) ^ 2;
t169 = xDDP(3);
t170 = xDDP(2);
t245 = t145 * t157;
t194 = (-t142 * t169 - t170 * t245) * t184;
t182 = 0.1e1 / pkin(3);
t251 = t182 * t282;
t37 = -(-t148 * t178 - t151 * t200) * pkin(2) + t55 * pkin(3);
t227 = t37 * t251;
t268 = t55 * t282;
t34 = (t227 - t268) * t184;
t239 = t34 * t37 * t282;
t209 = t185 * t239;
t203 = t142 * t148 + t145 * t151;
t25 = t34 * pkin(3) - t203 * t268;
t267 = t55 * t71;
t223 = t145 * t285;
t171 = xDDP(1);
t242 = t171 * t184;
t49 = t223 * t242;
t19 = -t25 * t185 * t267 + t49 + (t194 + t209) * t282;
t277 = pkin(2) * t19;
t244 = t146 * t158;
t193 = (-t143 * t169 - t170 * t244) * t184;
t250 = t182 * t281;
t38 = -(-t149 * t178 - t152 * t199) * pkin(2) + t56 * pkin(3);
t226 = t38 * t250;
t266 = t56 * t281;
t35 = (t226 - t266) * t184;
t238 = t35 * t38 * t281;
t208 = t185 * t238;
t202 = t143 * t149 + t146 * t152;
t26 = t35 * pkin(3) - t202 * t266;
t265 = t56 * t73;
t221 = t146 * t284;
t50 = t221 * t242;
t20 = -t26 * t185 * t265 + t50 + (t193 + t208) * t281;
t276 = pkin(2) * t20;
t243 = t147 * t159;
t192 = (-t144 * t169 - t170 * t243) * t184;
t249 = t182 * t280;
t39 = -(-t150 * t178 - t153 * t198) * pkin(2) + t57 * pkin(3);
t225 = t39 * t249;
t264 = t57 * t280;
t36 = (t225 - t264) * t184;
t237 = t36 * t39 * t280;
t207 = t185 * t237;
t263 = t57 * t75;
t201 = t144 * t150 + t147 * t153;
t27 = t36 * pkin(3) - t201 * t264;
t219 = t147 * t283;
t51 = t219 * t242;
t21 = -t27 * t185 * t263 + t51 + (t192 + t207) * t280;
t275 = pkin(2) * t21;
t137 = -t166 + t154;
t274 = sin(t137) / 0.2e1;
t139 = -t167 + t155;
t273 = sin(t139) / 0.2e1;
t141 = -t168 + t156;
t272 = sin(t141) / 0.2e1;
t136 = t166 + t154;
t271 = cos(t136) / 0.2e1;
t138 = t167 + t155;
t270 = cos(t138) / 0.2e1;
t140 = t168 + t156;
t269 = cos(t140) / 0.2e1;
t262 = t282 * (pkin(2) * t148 + pkin(3) * t142);
t261 = t281 * (pkin(2) * t149 + pkin(3) * t143);
t260 = t280 * (pkin(2) * t150 + pkin(3) * t144);
t259 = t170 - g(2);
t258 = t171 - g(1);
t257 = t142 * t282;
t256 = t143 * t281;
t255 = t144 * t280;
t181 = pkin(3) ^ 2;
t240 = 0.2e1 * pkin(3);
t31 = (-t268 + t227 / 0.2e1) * t184;
t254 = t182 * (-t34 * t181 + (-t203 * t31 * t240 + t268) * pkin(2));
t32 = (-t266 + t226 / 0.2e1) * t184;
t253 = t182 * (-t35 * t181 + (-t202 * t32 * t240 + t266) * pkin(2));
t33 = (-t264 + t225 / 0.2e1) * t184;
t252 = t182 * (-t36 * t181 + (-t201 * t33 * t240 + t264) * pkin(2));
t241 = t182 * t184;
t88 = pkin(2) * t151 + pkin(3) * t145;
t236 = t157 * t282 * t88;
t89 = pkin(2) * t152 + pkin(3) * t146;
t235 = t158 * t281 * t89;
t90 = pkin(2) * t153 + pkin(3) * t147;
t234 = t159 * t280 * t90;
t233 = t88 * t285;
t232 = t89 * t284;
t231 = t90 * t283;
t230 = t184 * t55 ^ 2 * t71;
t229 = t184 * t56 ^ 2 * t73;
t228 = t184 * t57 ^ 2 * t75;
t224 = t282 * t245;
t222 = t281 * t244;
t220 = t280 * t243;
t218 = t169 * t241;
t217 = t170 * t241;
t216 = t171 * t241;
t215 = t282 * t239;
t214 = t281 * t238;
t213 = t280 * t237;
t212 = -(t203 * pkin(2) + pkin(3)) * t209 * t251 - t216 * t233 + t217 * t236 + t218 * t262;
t211 = -(t202 * pkin(2) + pkin(3)) * t208 * t250 - t216 * t232 + t217 * t235 + t218 * t261;
t210 = -(t201 * pkin(2) + pkin(3)) * t207 * t249 - t216 * t231 + t217 * t234 + t218 * t260;
t121 = sin(t136);
t128 = cos(t137);
t206 = g(1) * t274 + g(2) * t271 + g(3) * t145 + t121 * t279 + t128 * t278;
t123 = sin(t138);
t130 = cos(t139);
t205 = g(1) * t273 + g(2) * t270 + g(3) * t146 + t123 * t279 + t130 * t278;
t125 = sin(t140);
t132 = cos(t141);
t204 = g(1) * t272 + g(2) * t269 + g(3) * t147 + t125 * t279 + t132 * t278;
t197 = g(1) * t271 + g(2) * t274 - g(3) * t142 + t121 * t278 + t128 * t279;
t196 = g(1) * t270 + g(2) * t273 - g(3) * t143 + t123 * t278 + t130 * t279;
t195 = g(1) * t269 + g(2) * t272 - g(3) * t144 + t125 * t278 + t132 * t279;
t191 = t282 * t194;
t190 = t281 * t193;
t189 = t280 * t192;
t183 = 0.1e1 / pkin(3) ^ 2;
t177 = cos(qJ(3,1));
t176 = cos(qJ(3,2));
t175 = cos(qJ(3,3));
t174 = sin(qJ(3,1));
t173 = sin(qJ(3,2));
t172 = sin(qJ(3,3));
t93 = t162 * g(1) - t159 * g(2);
t92 = t161 * g(1) - t158 * g(2);
t91 = t160 * g(1) - t157 * g(2);
t66 = -g(3) * t150 + t93 * t153;
t65 = g(3) * t153 + t93 * t150;
t64 = -g(3) * t149 + t92 * t152;
t63 = g(3) * t152 + t92 * t149;
t62 = -g(3) * t148 + t91 * t151;
t61 = g(3) * t151 + t91 * t148;
t60 = t258 * t159 + t259 * t162;
t59 = t258 * t158 + t259 * t161;
t58 = t258 * t157 + t259 * t160;
t18 = -t174 * t275 + t177 * t228 + t195;
t17 = -t173 * t276 + t176 * t229 + t196;
t16 = -t172 * t277 + t175 * t230 + t197;
t15 = t174 * t228 + t177 * t275 + t204;
t14 = t173 * t229 + t176 * t276 + t205;
t13 = t172 * t230 + t175 * t277 + t206;
t12 = t51 + t189 + (t213 + (-t27 - t252) * t263) * t185 + t210;
t11 = 0.2e1 * t51 + 0.2e1 * t189 + (0.2e1 * t213 + (-0.2e1 * t27 - t252) * t263) * t185 + t210;
t10 = t50 + t190 + (t214 + (-t26 - t253) * t265) * t185 + t211;
t9 = 0.2e1 * t50 + 0.2e1 * t190 + (0.2e1 * t214 + (-0.2e1 * t26 - t253) * t265) * t185 + t211;
t8 = t49 + t191 + (t215 + (-t25 - t254) * t267) * t185 + t212;
t7 = 0.2e1 * t49 + 0.2e1 * t191 + (0.2e1 * t215 + (-0.2e1 * t25 - t254) * t267) * t185 + t212;
t6 = -(t174 * t11 + (t183 * t75 * t39 - 0.2e1 * t249 * t264) * t39 * t177 * t185) * pkin(2) + t195;
t5 = -(t173 * t9 + (t183 * t73 * t38 - 0.2e1 * t250 * t266) * t38 * t176 * t185) * pkin(2) + t196;
t4 = -(t172 * t7 + (t183 * t71 * t37 - 0.2e1 * t251 * t268) * t37 * t175 * t185) * pkin(2) + t197;
t3 = pkin(2) * (t174 * t33 * t225 * t286 + t11 * t177) + t204;
t2 = pkin(2) * (t173 * t32 * t226 * t286 + t9 * t176) + t205;
t1 = pkin(2) * (t172 * t31 * t227 * t286 + t7 * t175) + t206;
t22 = [t157 * t58 + t158 * t59 + t159 * t60, (t19 * t223 + t20 * t221 + t21 * t219) * t184, (t65 * t219 + t63 * t221 + t61 * t223) * t184, (t66 * t219 + t64 * t221 + t62 * t223) * t184, (t10 * t221 + t12 * t219 + t8 * t223 + (-t10 * t232 - t12 * t231 - t8 * t233) * t182) * t184, (t1 * t223 + t2 * t221 + t3 * t219 + (-t13 * t233 - t14 * t232 - t15 * t231) * t182) * t184, (t4 * t223 + t5 * t221 + t6 * t219 + (-t16 * t233 - t17 * t232 - t18 * t231) * t182) * t184, t258; t160 * t58 + t161 * t59 + t162 * t60, (-t19 * t224 - t20 * t222 - t21 * t220) * t184, (-t65 * t220 - t63 * t222 - t61 * t224) * t184, (-t66 * t220 - t64 * t222 - t62 * t224) * t184, (-t10 * t222 - t12 * t220 - t8 * t224 + (t10 * t235 + t12 * t234 + t8 * t236) * t182) * t184, (-t1 * t224 - t2 * t222 - t3 * t220 + (t13 * t236 + t14 * t235 + t15 * t234) * t182) * t184, (-t4 * t224 - t5 * t222 - t6 * t220 + (t16 * t236 + t17 * t235 + t18 * t234) * t182) * t184, t259; 0, (-t19 * t257 - t20 * t256 - t21 * t255) * t184, (-t65 * t255 - t63 * t256 - t61 * t257) * t184, (-t66 * t255 - t64 * t256 - t62 * t257) * t184, (-t10 * t256 - t12 * t255 - t8 * t257 + (t10 * t261 + t12 * t260 + t8 * t262) * t182) * t184, (-t1 * t257 - t2 * t256 - t3 * t255 + (t13 * t262 + t14 * t261 + t15 * t260) * t182) * t184, (-t4 * t257 - t5 * t256 - t6 * t255 + (t16 * t262 + t17 * t261 + t18 * t260) * t182) * t184, t169 - g(3);];
tauX_reg  = t22;
