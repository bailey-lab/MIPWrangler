#pragma once
/*
 * mipUtils.hpp
 *
 *  Created on: Jul 26, 2016
 *      Author: nick
 */


#include "mipster/common.h"
#include "mipster/mipUtils/MipNameSorter.hpp"

namespace njhseq {


/**@brief get barcode number and read number, assumes name ends with _RX_BX_
 *
 * @param name name to process
 * @param barNum a reference to set to barcode number
 * @param readNum a reference to set to read number
 */
void processNameForBarReadCounts(const std::string & name,
		uint32_t & barNum, uint32_t & readNum);

/**@brief get read number, assumes name ends with _RX_BX_
 *
 * @param name name to process
 * @return the read count
 */
uint32_t getReadCntFromMipName(const std::string & name);

/**@brief get barcode number, assumes name ends with _RX_BX_
 *
 * @param name name to process
 * @return the barcode count
 */
uint32_t getBarcodeCntFromMipName(const std::string & name);


}  // namespace njhseq


