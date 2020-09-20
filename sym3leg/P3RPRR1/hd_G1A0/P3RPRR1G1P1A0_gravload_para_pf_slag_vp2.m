% Calculate Gravitation load for parallel robot
% P3RPRR1G1P1A0
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
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
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
% mrSges [4x3]
%   first moment of all robot links (mass times center of mass in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
%
% Output:
% taugX [3x1]
%   forces required to compensate gravitation load
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:23
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RPRR1G1P1A0_gravload_para_pf_slag_vp2(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRR1G1P1A0_gravload_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRR1G1P1A0_gravload_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRR1G1P1A0_gravload_para_pf_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRR1G1P1A0_gravload_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRR1G1P1A0_gravload_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RPRR1G1P1A0_gravload_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRR1G1P1A0_gravload_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRR1G1P1A0_gravload_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:23:15
% EndTime: 2020-03-09 21:23:15
% DurationCPUTime: 0.19s
% Computational Cost: add. (380->89), mult. (431->117), div. (18->4), fcn. (288->48), ass. (0->55)
t425 = legFrame(3,3);
t415 = sin(t425);
t418 = cos(t425);
t384 = -t415 * g(1) + t418 * g(2);
t387 = t418 * g(1) + t415 * g(2);
t429 = pkin(7) + qJ(3,3);
t412 = qJ(1,3) + t429;
t375 = -cos(t412) * (mrSges(3,1) * t384 - mrSges(3,2) * t387) + sin(t412) * (mrSges(3,1) * t387 + mrSges(3,2) * t384);
t390 = 0.1e1 / (pkin(1) * sin(t429) + sin(qJ(3,3)) * pkin(2));
t437 = t375 * t390;
t426 = legFrame(2,3);
t416 = sin(t426);
t419 = cos(t426);
t385 = -t416 * g(1) + t419 * g(2);
t388 = t419 * g(1) + t416 * g(2);
t430 = pkin(7) + qJ(3,2);
t413 = qJ(1,2) + t430;
t376 = -cos(t413) * (mrSges(3,1) * t385 - mrSges(3,2) * t388) + sin(t413) * (mrSges(3,1) * t388 + mrSges(3,2) * t385);
t391 = 0.1e1 / (pkin(1) * sin(t430) + sin(qJ(3,2)) * pkin(2));
t436 = t376 * t391;
t427 = legFrame(1,3);
t417 = sin(t427);
t420 = cos(t427);
t386 = -t417 * g(1) + t420 * g(2);
t389 = t420 * g(1) + t417 * g(2);
t431 = pkin(7) + qJ(3,1);
t414 = qJ(1,1) + t431;
t377 = -cos(t414) * (mrSges(3,1) * t386 - mrSges(3,2) * t389) + sin(t414) * (mrSges(3,1) * t389 + mrSges(3,2) * t386);
t392 = 0.1e1 / (pkin(1) * sin(t431) + sin(qJ(3,1)) * pkin(2));
t435 = t377 * t392;
t402 = mrSges(1,1) + (m(2) + m(3)) * pkin(1);
t421 = m(3) * pkin(2) + mrSges(2,1);
t422 = qJ(1,3) + pkin(7);
t434 = t390 * ((t387 * mrSges(2,2) - t384 * t421) * cos(t422) + (t384 * mrSges(2,2) + t421 * t387) * sin(t422) + (mrSges(1,2) * t387 - t384 * t402) * cos(qJ(1,3)) + (t384 * mrSges(1,2) + t402 * t387) * sin(qJ(1,3)) + t375);
t423 = qJ(1,2) + pkin(7);
t433 = t391 * ((t388 * mrSges(2,2) - t385 * t421) * cos(t423) + (t385 * mrSges(2,2) + t421 * t388) * sin(t423) + (mrSges(1,2) * t388 - t385 * t402) * cos(qJ(1,2)) + (t385 * mrSges(1,2) + t402 * t388) * sin(qJ(1,2)) + t376);
t424 = qJ(1,1) + pkin(7);
t432 = t392 * ((t389 * mrSges(2,2) - t386 * t421) * cos(t424) + (t386 * mrSges(2,2) + t421 * t389) * sin(t424) + (mrSges(1,2) * t389 - t386 * t402) * cos(qJ(1,1)) + (t386 * mrSges(1,2) + t402 * t389) * sin(qJ(1,1)) + t377);
t405 = t427 + t424;
t403 = t425 + t422;
t428 = 0.1e1 / pkin(3);
t411 = t427 + qJ(1,1);
t410 = t426 + qJ(1,2);
t409 = t425 + qJ(1,3);
t404 = t426 + t423;
t401 = qJ(3,1) + t405;
t400 = t426 + t413;
t399 = qJ(3,3) + t403;
t398 = cos(t401);
t397 = cos(t400);
t396 = cos(t399);
t395 = sin(t401);
t394 = sin(t400);
t393 = sin(t399);
t1 = [t396 * t434 + t397 * t433 + t398 * t432 - g(1) * m(4) + ((-pkin(1) * cos(t411) - pkin(2) * cos(t405) - pkin(3) * t398) * t435 + (-pkin(1) * cos(t410) - pkin(2) * cos(t404) - pkin(3) * t397) * t436 + (-pkin(1) * cos(t409) - pkin(2) * cos(t403) - pkin(3) * t396) * t437) * t428; t393 * t434 + t394 * t433 + t395 * t432 - g(2) * m(4) + ((-pkin(1) * sin(t411) - pkin(2) * sin(t405) - pkin(3) * t395) * t435 + (-pkin(1) * sin(t410) - pkin(2) * sin(t404) - pkin(3) * t394) * t436 + (-pkin(1) * sin(t409) - pkin(2) * sin(t403) - pkin(3) * t393) * t437) * t428; (-0.3e1 * m(2) - 0.3e1 * m(3) - m(4)) * g(3);];
taugX  = t1;
