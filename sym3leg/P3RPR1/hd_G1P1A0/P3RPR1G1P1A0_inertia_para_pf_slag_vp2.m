% Calculate inertia matrix for parallel robot
% P3RPR1G1P1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% qJ [2x3]
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
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d3]';
% m [3x1]
%   mass of all robot links (leg links until cut joint, platform)
% mrSges [3x3]
%   first moment of all robot links (mass times center of mass in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Ifges [3x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% MX [3x3]
%   inertia matrix in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-03 14:58
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3RPR1G1P1A0_inertia_para_pf_slag_vp2(xP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(2,3),zeros(3,3),zeros(3,3),zeros(4,1),zeros(2+1,1),zeros(2+1,3),zeros(2+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPR1G1P1A0_inertia_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [2 3]), ...
  'P3RPR1G1P1A0_inertia_para_pf_slag_vp2: qJ has to be [2x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RPR1G1P1A0_inertia_para_pf_slag_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [3 1]), ...
  'P3RPR1G1P1A0_inertia_para_pf_slag_vp2: m has to be [3x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [3,3]), ...
  'P3RPR1G1P1A0_inertia_para_pf_slag_vp2: mrSges has to be [3x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [3 6]), ...
  'P3RPR1G1P1A0_inertia_para_pf_slag_vp2: Ifges has to be [3x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPR1G1P1A0_inertia_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPR1G1P1A0_inertia_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 14:58:07
% EndTime: 2019-05-03 14:58:07
% DurationCPUTime: 0.29s
% Computational Cost: add. (1109->115), mult. (1601->204), div. (216->3), fcn. (1304->14), ass. (0->88)
t246 = 2 * mrSges(2,3);
t245 = 2 * pkin(1) * mrSges(2,1) + Ifges(2,2) + Ifges(1,3);
t230 = xP(3);
t217 = sin(t230);
t218 = cos(t230);
t231 = mrSges(3,2);
t232 = mrSges(3,1);
t244 = -t217 * t231 + t218 * t232;
t243 = -t217 * t232 - t218 * t231;
t242 = pkin(1) ^ 2;
t241 = koppelP(1,1);
t240 = koppelP(2,1);
t239 = koppelP(3,1);
t238 = koppelP(1,2);
t237 = koppelP(2,2);
t236 = koppelP(3,2);
t235 = 0.1e1 / qJ(2,1);
t234 = 0.1e1 / qJ(2,2);
t233 = 0.1e1 / qJ(2,3);
t229 = pkin(1) + pkin(2);
t228 = cos(qJ(1,1));
t227 = cos(qJ(1,2));
t226 = cos(qJ(1,3));
t225 = sin(qJ(1,1));
t224 = sin(qJ(1,2));
t223 = sin(qJ(1,3));
t222 = legFrame(1,3);
t221 = legFrame(2,3);
t220 = legFrame(3,3);
t216 = -m(2) * pkin(1) - mrSges(2,1);
t215 = cos(t222);
t214 = cos(t221);
t213 = cos(t220);
t212 = sin(t222);
t211 = sin(t221);
t210 = sin(t220);
t208 = t225 * qJ(2,1) + t229 * t228;
t207 = t224 * qJ(2,2) + t229 * t227;
t206 = t223 * qJ(2,3) + t229 * t226;
t205 = -t228 * qJ(2,1) + t225 * t229;
t204 = -t227 * qJ(2,2) + t224 * t229;
t203 = -t226 * qJ(2,3) + t223 * t229;
t202 = -t217 * t238 + t218 * t241;
t201 = -t217 * t237 + t218 * t240;
t200 = -t217 * t236 + t218 * t239;
t199 = -t217 * t241 - t218 * t238;
t198 = -t217 * t240 - t218 * t237;
t197 = -t217 * t239 - t218 * t236;
t196 = m(2) * (qJ(2,1) ^ 2 + t242) + qJ(2,1) * t246 + t245;
t195 = m(2) * (qJ(2,2) ^ 2 + t242) + qJ(2,2) * t246 + t245;
t194 = m(2) * (qJ(2,3) ^ 2 + t242) + qJ(2,3) * t246 + t245;
t193 = -t212 * t225 + t215 * t228;
t192 = t212 * t228 + t215 * t225;
t191 = -t211 * t224 + t214 * t227;
t190 = t211 * t227 + t214 * t224;
t189 = -t210 * t223 + t213 * t226;
t188 = t210 * t226 + t213 * t223;
t187 = -t212 * t205 + t208 * t215;
t186 = -t211 * t204 + t207 * t214;
t185 = -t210 * t203 + t206 * t213;
t184 = t205 * t215 + t212 * t208;
t183 = t204 * t214 + t211 * t207;
t182 = t203 * t213 + t210 * t206;
t181 = (m(2) * t187 + t193 * t216) * t235;
t180 = (m(2) * t184 + t192 * t216) * t235;
t179 = (m(2) * t186 + t191 * t216) * t234;
t178 = (m(2) * t183 + t190 * t216) * t234;
t177 = (m(2) * t185 + t189 * t216) * t233;
t176 = (m(2) * t182 + t188 * t216) * t233;
t175 = (t192 * t202 + t193 * t199) * t235;
t174 = (t190 * t201 + t191 * t198) * t234;
t173 = (t188 * t200 + t189 * t197) * t233;
t172 = (t187 * t216 + t193 * t196) * t235;
t171 = (t184 * t216 + t192 * t196) * t235;
t170 = (t186 * t216 + t191 * t195) * t234;
t169 = (t183 * t216 + t190 * t195) * t234;
t168 = (t185 * t216 + t189 * t194) * t233;
t167 = (t182 * t216 + t188 * t194) * t233;
t166 = (t184 * t202 + t187 * t199) * t235;
t165 = (t183 * t201 + t186 * t198) * t234;
t164 = (t182 * t200 + t185 * t197) * t233;
t163 = t166 * m(2) + t175 * t216;
t162 = t165 * m(2) + t174 * t216;
t161 = t164 * m(2) + t173 * t216;
t160 = t166 * t216 + t175 * t196;
t159 = t165 * t216 + t174 * t195;
t158 = t164 * t216 + t173 * t194;
t1 = [m(3) + (t172 * t193 + t181 * t187) * t235 + (t170 * t191 + t179 * t186) * t234 + (t168 * t189 + t177 * t185) * t233, (t172 * t192 + t181 * t184) * t235 + (t170 * t190 + t179 * t183) * t234 + (t168 * t188 + t177 * t182) * t233, t177 * t164 + t179 * t165 + t181 * t166 + t168 * t173 + t170 * t174 + t172 * t175 + t243; (t171 * t193 + t180 * t187) * t235 + (t169 * t191 + t178 * t186) * t234 + (t167 * t189 + t176 * t185) * t233, m(3) + (t171 * t192 + t180 * t184) * t235 + (t169 * t190 + t178 * t183) * t234 + (t167 * t188 + t176 * t182) * t233, t176 * t164 + t178 * t165 + t180 * t166 + t167 * t173 + t169 * t174 + t171 * t175 + t244; (t160 * t193 + t163 * t187) * t235 + (t159 * t191 + t162 * t186) * t234 + (t158 * t189 + t161 * t185) * t233 + t243, (t160 * t192 + t163 * t184) * t235 + (t159 * t190 + t162 * t183) * t234 + (t158 * t188 + t161 * t182) * t233 + t244, t158 * t173 + t159 * t174 + t160 * t175 + t161 * t164 + t162 * t165 + t163 * t166 + Ifges(3,3);];
MX  = t1;
