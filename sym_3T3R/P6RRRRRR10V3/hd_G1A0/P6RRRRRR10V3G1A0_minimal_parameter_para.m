% Return the minimum parameter vector for
% P6RRRRRR10V3G1A0
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,alpha3,alpha4,d1,d2,d3,d4]';
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
% koppelP [6x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% 
% Output:
% MPV [27x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-09-19 09:21
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = P6RRRRRR10V3G1A0_minimal_parameter_para(pkin, m, mrSges, Ifges, koppelP)

%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6),zeros(6,3)}
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'P6RRRRRR10V3G1A0_minimal_parameter_para: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P6RRRRRR10V3G1A0_minimal_parameter_para: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P6RRRRRR10V3G1A0_minimal_parameter_para: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P6RRRRRR10V3G1A0_minimal_parameter_para: Ifges has to be [4x6] (double)'); 
assert(isreal(koppelP) && all(size(koppelP) == [6 3]), ...
  'P6RRRRRR10V3G1A0_minimal_parameter_para: Koppelpunkt has to be [6x3] (double)');

%% Symbolic Calculation
% From minimal_parameter_parrob_matlab.m
unknown=NaN(27,1);
t1 = sin(pkin(4));
t2 = t1 ^ 2;
t3 = cos(pkin(5));
t4 = t3 ^ 2;
t5 = Ifges(3,2) * t4;
t8 = 0.2e1 * mrSges(3,3) * t4 * pkin(9);
t9 = pkin(2) ^ 2;
t10 = pkin(9) ^ 2;
t13 = m(3) * (t4 * t10 + t9);
t18 = t3 * pkin(9);
t20 = m(3) * t18 + mrSges(3,3) * t3 + mrSges(2,3);
t23 = pkin(1) ^ 2;
t24 = pkin(8) ^ 2;
t27 = m(2) + m(3);
t40 = sin(pkin(5));
t43 = pkin(2) * pkin(9);
t44 = m(3) * t40;
t54 = mrSges(3,3) * t40;
t60 = t40 ^ 2;
unknown(1,1) = Ifges(1,3) + (Ifges(2,2) + t5 + t8 + t13) * t2 + 0.2e1 * t20 * t2 * pkin(8) + t27 * (t2 * t24 + t23);
unknown(2,1) = t27 * pkin(1) + mrSges(1,1);
unknown(3,1) = -t27 * t1 * pkin(8) - t20 * t1 + mrSges(1,2);
unknown(4,1) = m(3) * t10 + 0.2e1 * pkin(9) * mrSges(3,3) + Ifges(2,1) - Ifges(2,2) + Ifges(3,2) - t13 - t5 - t8;
unknown(5,1) = mrSges(3,3) * t40 * pkin(2) + t44 * t43 + Ifges(2,4);
unknown(6,1) = -m(3) * t3 * t43 - mrSges(3,3) * t3 * pkin(2) + Ifges(2,5);
unknown(7,1) = Ifges(3,2) * t40 * t3 + t44 * t3 * t10 + 0.2e1 * t54 * t18 + Ifges(2,6);
unknown(8,1) = Ifges(2,3) + Ifges(3,2) * t60 + 0.2e1 * mrSges(3,3) * t60 * pkin(9) + m(3) * (t60 * t10 + t9);
unknown(9,1) = m(3) * pkin(2) + mrSges(2,1);
unknown(10,1) = -m(3) * t40 * pkin(9) + mrSges(2,2) - t54;
unknown(11,1) = Ifges(3,1) - Ifges(3,2);
unknown(12,1) = Ifges(3,4);
unknown(13,1) = Ifges(3,5);
unknown(14,1) = Ifges(3,6);
unknown(15,1) = Ifges(3,3);
unknown(16,1) = mrSges(3,1);
unknown(17,1) = mrSges(3,2);
unknown(18,1) = Ifges(4,1);
unknown(19,1) = Ifges(4,4);
unknown(20,1) = Ifges(4,5);
unknown(21,1) = Ifges(4,2);
unknown(22,1) = Ifges(4,6);
unknown(23,1) = Ifges(4,3);
unknown(24,1) = mrSges(4,1);
unknown(25,1) = mrSges(4,2);
unknown(26,1) = mrSges(4,3);
unknown(27,1) = m(4);
MPV  = unknown;
