% Calculate inertia matrix for parallel robot
% P3PRRR2G3P1A0
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
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4]';
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
% MX [3x3]
%   inertia matrix in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:20
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3PRRR2G3P1A0_inertia_para_pf_slag_vp1(xP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR2G3P1A0_inertia_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR2G3P1A0_inertia_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRR2G3P1A0_inertia_para_pf_slag_vp1: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRR2G3P1A0_inertia_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3PRRR2G3P1A0_inertia_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3PRRR2G3P1A0_inertia_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR2G3P1A0_inertia_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR2G3P1A0_inertia_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:20:01
% EndTime: 2020-03-09 21:20:01
% DurationCPUTime: 0.24s
% Computational Cost: add. (481->65), mult. (579->123), div. (180->5), fcn. (252->18), ass. (0->71)
t250 = pkin(1) * m(3);
t230 = rSges(3,1) ^ 2 + rSges(3,2) ^ 2;
t194 = t230 * m(3) + Icges(3,3);
t216 = sin(qJ(3,3));
t228 = (rSges(3,1) * cos(qJ(3,3)) - rSges(3,2) * t216) * t250;
t185 = t228 + t194;
t224 = 0.1e1 / pkin(2);
t249 = t185 * t224;
t217 = sin(qJ(3,2));
t227 = (rSges(3,1) * cos(qJ(3,2)) - rSges(3,2) * t217) * t250;
t186 = t227 + t194;
t248 = t186 * t224;
t218 = sin(qJ(3,1));
t226 = (rSges(3,1) * cos(qJ(3,1)) - rSges(3,2) * t218) * t250;
t187 = t226 + t194;
t247 = t187 * t224;
t207 = -legFrame(3,2) + qJ(2,3);
t202 = qJ(3,3) + t207;
t195 = sin(t202);
t188 = pkin(2) * t195 + pkin(1) * sin(t207);
t210 = 0.1e1 / t216;
t246 = t188 * t210;
t208 = -legFrame(2,2) + qJ(2,2);
t203 = qJ(3,2) + t208;
t196 = sin(t203);
t189 = pkin(2) * t196 + pkin(1) * sin(t208);
t211 = 0.1e1 / t217;
t245 = t189 * t211;
t209 = -legFrame(1,2) + qJ(2,1);
t204 = qJ(3,1) + t209;
t197 = sin(t204);
t190 = pkin(2) * t197 + pkin(1) * sin(t209);
t212 = 0.1e1 / t218;
t244 = t190 * t212;
t198 = cos(t202);
t191 = -pkin(2) * t198 - pkin(1) * cos(t207);
t243 = t191 * t210;
t199 = cos(t203);
t192 = -pkin(2) * t199 - pkin(1) * cos(t208);
t242 = t192 * t211;
t200 = cos(t204);
t193 = -pkin(2) * t200 - pkin(1) * cos(t209);
t241 = t193 * t212;
t240 = t194 * t224;
t239 = t195 * t210;
t238 = t196 * t211;
t237 = t197 * t212;
t236 = t198 * t210;
t235 = t199 * t211;
t234 = t200 * t212;
t225 = 0.1e1 / pkin(1);
t233 = t210 * t225;
t232 = t211 * t225;
t231 = t212 * t225;
t229 = Icges(2,3) + Icges(3,3) + (pkin(1) ^ 2 + t230) * m(3) + (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2);
t184 = 0.2e1 * t226 + t229;
t183 = 0.2e1 * t227 + t229;
t182 = 0.2e1 * t228 + t229;
t181 = (t187 * t200 + t193 * t240) * t231;
t180 = (t186 * t199 + t192 * t240) * t232;
t179 = (t185 * t198 + t191 * t240) * t233;
t178 = (-t187 * t197 + t190 * t240) * t231;
t177 = (-t186 * t196 + t189 * t240) * t232;
t176 = (-t185 * t195 + t188 * t240) * t233;
t175 = (t184 * t200 + t193 * t247) * t231;
t174 = (t183 * t199 + t192 * t248) * t232;
t173 = (t182 * t198 + t191 * t249) * t233;
t172 = (-t184 * t197 + t190 * t247) * t231;
t171 = (-t183 * t196 + t189 * t248) * t232;
t170 = (-t182 * t195 + t188 * t249) * t233;
t1 = [m(4) + (-t170 * t239 - t171 * t238 - t172 * t237 + (t176 * t246 + t177 * t245 + t178 * t244) * t224) * t225, (t170 * t236 + t171 * t235 + t172 * t234 + (t176 * t243 + t177 * t242 + t178 * t241) * t224) * t225, 0; (-t173 * t239 - t174 * t238 - t175 * t237 + (t179 * t246 + t180 * t245 + t181 * t244) * t224) * t225, m(4) + (t173 * t236 + t174 * t235 + t175 * t234 + (t179 * t243 + t180 * t242 + t181 * t241) * t224) * t225, 0; 0, 0, (3 * m(1)) + 0.3e1 * m(2) + 0.3e1 * m(3) + m(4);];
MX  = t1;
