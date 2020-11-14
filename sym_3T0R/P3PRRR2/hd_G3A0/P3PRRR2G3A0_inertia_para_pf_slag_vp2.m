% Calculate inertia matrix for parallel robot
% P3PRRR2G3A0
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
% MX [3x3]
%   inertia matrix in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:20
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3PRRR2G3A0_inertia_para_pf_slag_vp2(xP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR2G3A0_inertia_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR2G3A0_inertia_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRR2G3A0_inertia_para_pf_slag_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRR2G3A0_inertia_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3PRRR2G3A0_inertia_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3PRRR2G3A0_inertia_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR2G3A0_inertia_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR2G3A0_inertia_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:20:04
% EndTime: 2020-03-09 21:20:04
% DurationCPUTime: 0.16s
% Computational Cost: add. (361->60), mult. (339->116), div. (180->5), fcn. (252->18), ass. (0->68)
t205 = 0.1e1 / pkin(2);
t229 = Ifges(3,3) * t205;
t189 = -legFrame(3,2) + qJ(2,3);
t186 = qJ(3,3) + t189;
t180 = sin(t186);
t168 = pkin(2) * t180 + pkin(1) * sin(t189);
t198 = sin(qJ(3,3));
t192 = 0.1e1 / t198;
t228 = t168 * t192;
t190 = -legFrame(2,2) + qJ(2,2);
t187 = qJ(3,2) + t190;
t181 = sin(t187);
t169 = pkin(2) * t181 + pkin(1) * sin(t190);
t199 = sin(qJ(3,2));
t193 = 0.1e1 / t199;
t227 = t169 * t193;
t191 = -legFrame(1,2) + qJ(2,1);
t188 = qJ(3,1) + t191;
t182 = sin(t188);
t170 = pkin(2) * t182 + pkin(1) * sin(t191);
t200 = sin(qJ(3,1));
t194 = 0.1e1 / t200;
t226 = t170 * t194;
t183 = cos(t186);
t171 = -pkin(2) * t183 - pkin(1) * cos(t189);
t225 = t171 * t192;
t184 = cos(t187);
t172 = -pkin(2) * t184 - pkin(1) * cos(t190);
t224 = t172 * t193;
t185 = cos(t188);
t173 = -pkin(2) * t185 - pkin(1) * cos(t191);
t223 = t173 * t194;
t209 = (mrSges(3,1) * cos(qJ(3,3)) - mrSges(3,2) * t198) * pkin(1);
t177 = Ifges(3,3) + t209;
t222 = t177 * t205;
t208 = (mrSges(3,1) * cos(qJ(3,2)) - mrSges(3,2) * t199) * pkin(1);
t178 = Ifges(3,3) + t208;
t221 = t178 * t205;
t207 = (mrSges(3,1) * cos(qJ(3,1)) - mrSges(3,2) * t200) * pkin(1);
t179 = Ifges(3,3) + t207;
t220 = t179 * t205;
t219 = t180 * t192;
t218 = t181 * t193;
t217 = t182 * t194;
t216 = t183 * t192;
t215 = t184 * t193;
t214 = t185 * t194;
t206 = 0.1e1 / pkin(1);
t213 = t192 * t206;
t212 = t193 * t206;
t211 = t194 * t206;
t210 = m(3) * pkin(1) ^ 2 + Ifges(2,3) + Ifges(3,3);
t176 = 0.2e1 * t207 + t210;
t175 = 0.2e1 * t208 + t210;
t174 = 0.2e1 * t209 + t210;
t167 = (t173 * t229 + t179 * t185) * t211;
t166 = (t172 * t229 + t178 * t184) * t212;
t165 = (t171 * t229 + t177 * t183) * t213;
t164 = (t170 * t229 - t179 * t182) * t211;
t163 = (t169 * t229 - t178 * t181) * t212;
t162 = (t168 * t229 - t177 * t180) * t213;
t161 = (t173 * t220 + t176 * t185) * t211;
t160 = (t172 * t221 + t175 * t184) * t212;
t159 = (t171 * t222 + t174 * t183) * t213;
t158 = (t170 * t220 - t176 * t182) * t211;
t157 = (t169 * t221 - t175 * t181) * t212;
t156 = (t168 * t222 - t174 * t180) * t213;
t1 = [m(4) + (-t156 * t219 - t157 * t218 - t158 * t217 + (t162 * t228 + t163 * t227 + t164 * t226) * t205) * t206, (t156 * t216 + t157 * t215 + t158 * t214 + (t162 * t225 + t163 * t224 + t164 * t223) * t205) * t206, 0; (-t159 * t219 - t160 * t218 - t161 * t217 + (t165 * t228 + t166 * t227 + t167 * t226) * t205) * t206, m(4) + (t159 * t216 + t160 * t215 + t161 * t214 + (t165 * t225 + t166 * t224 + t167 * t223) * t205) * t206, 0; 0, 0, (3 * m(1)) + (3 * m(2)) + 0.3e1 * m(3) + m(4);];
MX  = t1;
